#include <iostream>

#include "utils/zmq_utils.h"


void send_message(zmq::socket_t *socket, std::vector<zmq::message_t> &msg){
    for (auto it = msg.begin(); it != msg.end(); ++it) {
        if (it + 1 != msg.end()) {
            socket->send(*it, ZMQ_SNDMORE);
        } else {
            socket->send(*it, 0);
        }
    }
}

bool recv_message(zmq::socket_t *socket, std::vector<zmq::message_t> &msg, bool wait){
    if (!msg.empty()){
        std::cerr << "接收消息错误，list中有内容" << std::endl;
        exit(-1);
    }
    msg.emplace_back();
    if (wait) {
        socket->recv(&msg[0], 0);
    } else {
        if (!socket->recv(&msg[0], ZMQ_DONTWAIT)) {
            msg.clear();
            return false;
        }
    }
    bool has_more = msg[0].more();
    while (has_more) {
        msg.emplace_back();
        socket->recv(&msg.back(), 0);
        has_more = msg.back().more();
    }
    return true;
}

void delete_socket(zmq::socket_t *socket) {
    int i = 0;
    if (socket != NULL) {
        socket->setsockopt(ZMQ_LINGER, &i, sizeof(i));
        delete socket;
    }
}


