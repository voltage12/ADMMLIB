#ifndef COORDINATOR_H
#define COORDINATOR_H

#include <vector>
#include <zmq.hpp>
#include <map>
#include <mutex>
#include <condition_variable>

#include "calculator.h"
#include "utils/properties.h"
#include "optimizer/sparse_sample_set.h"


class Coordinator {
    friend void worker_thread_function(Coordinator *coordinator, int id);
public:
    enum Role {Master, Leader};

    Coordinator(const std::string  &hostfile_path, Properties &properties);

    void run();

    inline int get_worker_num() { return worker_num; }

    virtual ~Coordinator();

private:
    Role role;
    bool is_stop;
    int dim; /* 变量维度 */
    int leader_id;
    int starting_id; /* 全局id */
    int worker_num; /* 总进程数 */
    int thread_num; /* 本节点线程数 */
    int min_barrier, max_delay;
    int z_clock, max_iter_num;
    double rho;
    double comm_time;
    double ABSTOL, RELTOL;
    double *z_old, *z_new, *rho_x_plus_y;
    timeval start_time;
    zmq::context_t context;
    zmq::socket_t *left_socket, *right_socket;
    zmq::socket_t *master_worker_socket, *master_leader_socket;
    std::string my_hostname;
    std::string master_hostname;
    std::vector<std::string> host_list;
    std::map<int, double> cal_time;
    std::map<int, double> wait_time;
    std::map<int, double *> x_ptrs;
    std::map<int, double *> y_ptrs;
    std::map<int, double *> z_ptrs;
    std::map<int, double *> x_caches;
    std::map<int, double *> y_caches;
    std::map<int, Calculator *> calculators;
    std::map<int, bool> thread_is_stop;
    std::map<int, bool> can_calculate;
    std::map<int, std::mutex *> mutexs;
    std::map<int, std::condition_variable *> cond_vars;
    std::map<int, int> worker_delay;
    SparseSampleSet *test_sample_set;

    void parse_hostfile(const std::string &hostfile_path);

    void act_master();

    void act_leader();

    bool check_worker_delay();

    void terminal();

    bool update_z(size_t ready_worker_num, const int *ready_worker_list);

    void ring_allreduce_sum(double *buf, int length);

    void allreduce_sum(double *buf, int length);

    void predict();

    double cal_objective_value();

    void barrier();

    void awake_waiting_worker(size_t ready_worker_num, const int *ready_worker_list);
};


#endif //MCCADMM_COORDINATOR_H
