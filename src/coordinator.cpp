#include <iostream>
#include <thread>
#include <assert.h>
#include <cmath>
#include <unistd.h>
#include <algorithm>
#include <sys/time.h>
#include <fstream>

#include "coordinator.h"
#include "lr_calculator.h"
#include "utils/zmq_utils.h"


void worker_thread_function(Coordinator *coordinator, int id) {
    zmq::socket_t *socket = new zmq::socket_t(coordinator->context, ZMQ_DEALER);
    std::string address = std::string("tcp://") + coordinator->master_hostname + ":50001";
    socket->connect(address);
    double rho = coordinator->rho;
    double *x = coordinator->x_ptrs[id];
    double *y = coordinator->y_ptrs[id];
    double *z = coordinator->z_ptrs[id];
    Calculator *calculator = coordinator->calculators[id];
    timeval start_time, end_time;
    timeval cal_start_time, cal_end_time;
    gettimeofday(&start_time, NULL);
    while (!coordinator->is_stop) {
        std::unique_lock<std::mutex> lck(*coordinator->mutexs[id]);
        while (!coordinator->can_calculate[id]) {
            coordinator->cond_vars[id]->wait(lck);
        }
        /* 计算线程被唤醒可能是因为算法已经结束，因此需要先检查条件 */
        if (coordinator->is_stop) {
            break;
        }
        gettimeofday(&cal_start_time, NULL);
        calculator->update_x(rho, y, z, x);
        calculator->update_y(rho, x, z, y);
        coordinator->can_calculate[id] = false;
        gettimeofday(&cal_end_time, NULL);
        coordinator->cal_time[id] += ((cal_end_time.tv_sec - cal_start_time.tv_sec) +
                                      (cal_end_time.tv_usec - cal_start_time.tv_usec) / 1000000.0);
        zmq::message_t msg(&id, sizeof(id));
        socket->send(msg, 0);
    }
    gettimeofday(&end_time, NULL);
    coordinator->wait_time[id] =
            (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0 -
            coordinator->cal_time[id];
    delete_socket(socket);
    coordinator->thread_is_stop[id] = true;
}

Coordinator::Coordinator(const std::string &hostfile_path, Properties &properties) {
    parse_hostfile(hostfile_path);
    is_stop = false;
    comm_time = 0;
    z_clock = 0;
    dim = properties.get_int("dim");
    min_barrier = properties.get_int("min_barrier");
    max_delay = properties.get_int("max_delay");
    max_iter_num = properties.get_int("max_iter_num");
    rho = properties.get_double("rho");
    ABSTOL = properties.get_double("ABSTOL");
    RELTOL = properties.get_double("RELTOL");
    for (int i = starting_id; i < starting_id + thread_num; ++i) {
        x_ptrs[i] = new double[dim];
        y_ptrs[i] = new double[dim];
        z_ptrs[i] = new double[dim];
        x_caches[i] = new double[dim];
        y_caches[i] = new double[dim];
        thread_is_stop[i] = false;
        can_calculate[i] = true;
        calculators[i] = new LRCalculator(i, properties);
        calculators[i]->init(x_ptrs[i], y_ptrs[i], z_ptrs[i]);
        memcpy(x_caches[i], x_ptrs[i], dim * sizeof(double));
        memcpy(y_caches[i], y_ptrs[i], dim * sizeof(double));
        cal_time[i] = 0;
        wait_time[i] = 0;
        mutexs[i] = new std::mutex();
        cond_vars[i] = new std::condition_variable();
    }
    z_old = new double[dim];
    z_new = new double[dim];
    rho_x_plus_y = new double[dim];
    memcpy(z_new, z_ptrs[starting_id], dim * sizeof(double));
    test_sample_set = NULL;
    if (role == Role::Master) {
        std::string test_file_path = properties.get_string("test_file_path");
        test_sample_set = new SparseSampleSet(test_file_path);
        properties.print();
    }
}

Coordinator::~Coordinator() {
    delete_socket(master_leader_socket);
    delete_socket(master_worker_socket);
    delete_socket(left_socket);
    delete_socket(right_socket);
    delete test_sample_set;
    delete[] z_old;
    delete[] z_new;
    delete[] rho_x_plus_y;
    for (int i = starting_id; i < starting_id + thread_num; ++i) {
        delete[] x_ptrs[i];
        delete[] y_ptrs[i];
        delete[] z_ptrs[i];
        delete[] x_caches[i];
        delete[] y_caches[i];
        delete calculators[i];
        delete mutexs[i];
        delete cond_vars[i];
    }
}

void Coordinator::parse_hostfile(const std::string &hostfile_path) {
    char temp[10];
    gethostname(temp, 10);
    my_hostname.assign(temp);
    std::ifstream reader(hostfile_path);
    if (reader.fail()) {
        std::cerr << "无法打开文件:" << hostfile_path << std::endl;
        exit(-1);
    }
    std::string hostname;
    std::map<std::string, int> temp_map;
    while (std::getline(reader, hostname)) {
        if (temp_map.count(hostname) == 0) {
            host_list.push_back(hostname);
            temp_map[hostname] = 1;
        } else {
            temp_map[hostname] += 1;
        }
    }
    reader.close();
    master_hostname = host_list.front();
    if (my_hostname == master_hostname) {
        role = Role::Master;
        master_leader_socket = new zmq::socket_t(context, ZMQ_ROUTER);
        master_leader_socket->bind("tcp://*:50000");
        master_worker_socket = new zmq::socket_t(context, ZMQ_ROUTER);
        master_worker_socket->bind("tcp://*:50001");
    } else {
        role = Role::Leader;
        master_leader_socket = new zmq::socket_t(context, ZMQ_DEALER);
        master_leader_socket->setsockopt(ZMQ_IDENTITY, my_hostname.data(), my_hostname.size());
        std::string address = std::string("tcp://") + master_hostname + ":50000";
        master_leader_socket->connect(address);
        master_worker_socket = NULL;
    }
    left_socket = new zmq::socket_t(context, ZMQ_ROUTER);
    left_socket->bind("tcp://*:50002");
    right_socket = new zmq::socket_t(context, ZMQ_DEALER);
    worker_num = 0;
    for (int i = 0; i < host_list.size(); ++i) {
        if (host_list[i] == my_hostname) {
            leader_id = i;
            starting_id = worker_num;
            thread_num = temp_map[host_list[i]];
            std::string address = std::string("tcp://") + host_list[(i + 1) % host_list.size()] + ":50002";
            right_socket->connect(address);
        }
        worker_num += temp_map[host_list[i]];
    }
    std::cout << worker_num << ' ' << leader_id << ' ' << starting_id << ' ' << thread_num << " My hostname is "
              << my_hostname
              << ", Master hostname is " << master_hostname << std::endl;
}

void Coordinator::run() {
    barrier();
    for (int i = starting_id; i < starting_id + thread_num; ++i) {
        std::thread t(worker_thread_function, this, i);
        t.detach();
    }
    if (role == Role::Master) {
        act_master();
    } else if (role == Role::Leader) {
        act_leader();
    }
    terminal();
    double send_buf[] = {0, 0};
    for (int i = starting_id; i < starting_id + thread_num; ++i) {
        send_buf[0] += cal_time[i];
        send_buf[1] += wait_time[i];
    }
    allreduce_sum(send_buf, 2);
    if (role == Role::Master) {
        std::cout << "总通信时间：" << comm_time << std::endl;
        std::cout << "平均计算时间：" << send_buf[0] / worker_num << "，平均等待时间：" << send_buf[1] / worker_num << std::endl;
        predict();
    }
}

void Coordinator::act_leader() {
    std::vector<zmq::message_t> msg;
    while (!is_stop) {
        recv_message(master_leader_socket, msg, true);
        if (update_z(msg.front().size() / sizeof(int), msg.front().data<int>())) {
            return;
        }
        awake_waiting_worker(msg.front().size() / sizeof(int), msg.front().data<int>());
        msg.clear();
    }
}

void Coordinator::act_master() {
    printf("%3s %10s %10s %10s %10s %10s %10s\n", "#", "r norm", "esp_pri", "s norm", "esp_dual", "obj_val", "time");
    fflush(stdout);
    std::vector<int> ready_worker_list;
    std::vector<zmq::message_t> msg;
    for (int i = 0; i < worker_num; ++i) {
        worker_delay[i] = 0;
    }
    gettimeofday(&start_time, NULL);
    while (!is_stop) {
        recv_message(master_worker_socket, msg, true);
        int worker_id = *msg[1].data<int>();
        worker_delay[worker_id] = -1;
        ready_worker_list.push_back(worker_id);
        msg.clear();
        if (ready_worker_list.size() >= min_barrier && check_worker_delay()) {
            std::sort(ready_worker_list.begin(), ready_worker_list.end());
            for (auto it = host_list.begin() + 1; it != host_list.end(); ++it) {
                msg.emplace_back(it->data(), it->size());
                msg.emplace_back(ready_worker_list.data(), ready_worker_list.size() * sizeof(int));
                send_message(master_leader_socket, msg);
                msg.clear();
            }
            if (update_z(ready_worker_list.size(), ready_worker_list.data())) {
                return;
            }
            awake_waiting_worker(ready_worker_list.size(), ready_worker_list.data());
            ready_worker_list.clear();
            for (auto it = worker_delay.begin(); it != worker_delay.end(); ++it) {
                ++it->second;
            }
        }
    }
}

bool Coordinator::update_z(size_t ready_worker_num, const int *ready_worker_list) {
    /* 首先将当前z变量备份，并且将rho.x + y初始化 */
    for (int i = 0; i < dim; ++i) {
        z_old[i] = z_new[i];
        rho_x_plus_y[i] = 0;
    }
    int a = starting_id, b = 0;
    while (a < starting_id + thread_num && b < ready_worker_num) {
        if (a == ready_worker_list[b]) {
            memcpy(x_caches[a], x_ptrs[a], dim * sizeof(double));
            memcpy(y_caches[a], y_ptrs[a], dim * sizeof(double));
            ++a;
            ++b;
            continue;
        }
        a > ready_worker_list[b] ? ++b : ++a;
    }
    /* 接着leader遍历组内各成员（包括自己）的缓存，计算组内的rho.x + y的累加和 */
    for (a = starting_id; a < starting_id + thread_num; ++a) {
        double *x_temp = x_caches[a];
        double *y_temp = y_caches[a];
        for (int i = 0; i < dim; ++i) {
            rho_x_plus_y[i] += (rho * x_temp[i] + y_temp[i]);
        }
    }
    ring_allreduce_sum(rho_x_plus_y, dim);
    calculators[starting_id]->update_z(worker_num, rho, rho_x_plus_y, z_new);
    ++z_clock;
    double reduce_buf[] = {0, 0, 0};
    /* 计算并累加组内成员的||x_i||_2^2), ||y_i||_2^2), ||r_i||_2^2) */
    for (a = starting_id; a < starting_id + thread_num; ++a) {
        double *x_temp = x_caches[a];
        double *y_temp = y_caches[a];
        for (int j = 0; j < dim; ++j) {
            reduce_buf[0] += x_temp[j] * x_temp[j];
            reduce_buf[1] += y_temp[j] * y_temp[j];
            double temp = x_temp[j] - z_new[j];
            reduce_buf[2] += temp * temp;
        }
    }
    allreduce_sum(reduce_buf, 3);
    double nxstack = sqrt(reduce_buf[0]); /* sqrt(sum ||x_i||_2^2) */
    double nystack = sqrt(reduce_buf[1]); /* sqrt(sum ||y_i||_2^2) */
    double prires = sqrt(reduce_buf[2]); /* sqrt(sum ||r_i||_2^2) */
    double z_diff = 0; /* 存放||z_new - z_old||_2^2 */
    double z_norm = 0; /* 存放||z_new||_2^2 */
    for (int i = 0; i < dim; ++i) {
        double temp = z_old[i] - z_new[i];
        z_diff += temp * temp;
        z_norm += z_new[i] * z_new[i];
    }
    double dualres = rho * sqrt(worker_num * z_diff);
    double eps_pri = sqrt(dim * worker_num) * ABSTOL + RELTOL * fmax(nxstack, sqrt(worker_num * z_norm));
    double eps_dual = sqrt(dim * worker_num) * ABSTOL + RELTOL * nystack;
    if (role == Role::Master) {
        timeval end_time;
        gettimeofday(&end_time, NULL);
        double wait_time =
                (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
        printf("%3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", z_clock, prires, eps_pri, dualres, eps_dual,
               cal_objective_value(), wait_time);
        fflush(stdout);
    }
    return ((prires <= eps_pri && dualres <= eps_dual) || z_clock >= max_iter_num);
}

void Coordinator::awake_waiting_worker(size_t ready_worker_num, const int *ready_worker_list) {
    int a = starting_id;
    int b = 0;
    while (a < starting_id + thread_num && b < ready_worker_num) {
        if (a == ready_worker_list[b]) {
            memcpy(z_ptrs[a], z_new, dim * sizeof(double));
            can_calculate[a] = true;
            cond_vars[a]->notify_one();
            ++a;
            ++b;
            continue;
        }
        a > ready_worker_list[b] ? ++b : ++a;
    }
}

bool Coordinator::check_worker_delay() {
    for (auto it = worker_delay.begin(); it != worker_delay.end(); ++it) {
        if (it->second >= max_delay) {
            std::cout << "worker" << it->first << "达到最大延迟，将等待其更新到达" << std::endl;
            return false;
        }
    }
    return true;
}

void Coordinator::ring_allreduce_sum(double *buf, int length) {
    std::vector<zmq::message_t> msg;
    int leader_num = static_cast<int>(host_list.size());
    if (leader_num == 1) {
        return;
    }
    if (length < leader_num * (leader_num - 1)) {
        std::cerr << "length is too small" << std::endl;
        exit(-1);
    }
    /* 待reduce数据被切分成leader_num个段 */
    int temp = static_cast<int>(std::ceil(1.0 * length / leader_num));
    /* seg编号 -> (seg起始地址，seg长度) */
    std::map<int, std::pair<int, int>> seg;
    for (int i = 0; i < leader_num; ++i) {
        seg[i].first = i * temp;
        seg[i].second = temp;
    }
    seg[leader_num - 1].second = (length - (leader_num - 1) * temp);
    int send_seg_index = (leader_id - 1 + leader_num) % leader_num;
    int recv_seg_index = (send_seg_index - 1 + leader_num) % leader_num;
    timeval comm_start_time;
    gettimeofday(&comm_start_time, NULL);
    for (int i = 0; i < leader_num - 1; ++i) {
        msg.emplace_back(buf + seg[send_seg_index].first, seg[send_seg_index].second * sizeof(double));
        send_message(right_socket, msg);
        msg.clear();
        recv_message(left_socket, msg, true);
        double *recv_buf = msg[1].data<double>();
        assert(msg[1].size() == seg[recv_seg_index].second * sizeof(double));
        double *base_ptr = buf + seg[recv_seg_index].first;
        double seg_length = seg[recv_seg_index].second;
        for (int j = 0; j < seg_length; ++j) {
            base_ptr[j] += recv_buf[j];
        }
        msg.clear();
        send_seg_index = recv_seg_index;
        recv_seg_index = (recv_seg_index - 1 + leader_num) % leader_num;
    }
    for (int i = 0; i < leader_num - 1; ++i) {
        msg.emplace_back(buf + seg[send_seg_index].first, seg[send_seg_index].second * sizeof(double));
        send_message(right_socket, msg);
        msg.clear();
        recv_message(left_socket, msg, true);
        assert(msg[1].size() == seg[recv_seg_index].second * sizeof(double));
        memcpy(buf + seg[recv_seg_index].first, msg[1].data(), msg[1].size());
        msg.clear();
        send_seg_index = recv_seg_index;
        recv_seg_index = (recv_seg_index - 1 + leader_num) % leader_num;
    }
    timeval end_time;
    gettimeofday(&end_time, NULL);
    comm_time += ((end_time.tv_sec - comm_start_time.tv_sec) +
                  (end_time.tv_usec - comm_start_time.tv_usec) / 1000000.0);
}

void Coordinator::allreduce_sum(double *buf, int length) {
    std::vector<zmq::message_t> msg;
    if (role == Role::Master) {
        int count = 1;
        while (count < host_list.size()) {
            recv_message(master_leader_socket, msg, true);
            double *temp = msg[1].data<double>();
            for (int i = 0; i < length; ++i) {
                buf[i] += temp[i];
            }
            msg.clear();
            ++count;
        }
        for (auto it = host_list.begin() + 1; it != host_list.end(); ++it) {
            msg.emplace_back(it->data(), it->size());
            msg.emplace_back(buf, length * sizeof(double));
            send_message(master_leader_socket, msg);
            msg.clear();
        }
    } else {
        msg.emplace_back(buf, length * sizeof(double));
        send_message(master_leader_socket, msg);
        msg.clear();
        recv_message(master_leader_socket, msg, true);
        memcpy(buf, msg.front().data(), msg.front().size());
        msg.clear();
    }
}

void Coordinator::predict() {
    int counter = 0;
    for (int i = 0; i < test_sample_set->get_sample_num(); ++i) {
        double temp = 1.0 / (1 + exp(-1 * test_sample_set->dot(i, z_new)));
        if (test_sample_set->get_label(i) == 1 && temp >= 0.5) {
            ++counter;
        }
        if (test_sample_set->get_label(i) == -1 && temp < 0.5) {
            ++counter;
        }
    }
    std::cout << "预测精度:" << counter * 100.0 / test_sample_set->get_sample_num() << "%" << std::endl;
}

double Coordinator::cal_objective_value() {
    double sum = 0;
    int sample_num = test_sample_set->get_sample_num();
    for (int i = 0; i < sample_num; ++i) {
        sum += std::log(1 + std::exp(-test_sample_set->get_label(i) * test_sample_set->dot(i, z_new)));
    }
    return sum;
}

void Coordinator::terminal() {
    is_stop = true;
    for (int i = starting_id; i < starting_id + thread_num; ++i) {
        can_calculate[i] = true;
        cond_vars[i]->notify_one();
    }
    bool flag;
    do {
        flag = true;
        for (int i = starting_id; i < starting_id + thread_num; ++i) {
            if (!thread_is_stop[i]) {
                flag = false;
                sleep(1);
                continue;
            }
        }
    } while (!flag);
}

void Coordinator::barrier() {
    double tag = -1;
    /* 借助allreduce实现barrier */
    allreduce_sum(&tag, 1);
}
