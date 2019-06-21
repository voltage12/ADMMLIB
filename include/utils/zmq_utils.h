#ifndef ZMQ_UTILS_H
#define ZMQ_UTILS_H

#include <zmq.hpp>
#include <vector>


void send_message(zmq::socket_t *socket, std::vector<zmq::message_t> &msg);

bool recv_message(zmq::socket_t *socket, std::vector<zmq::message_t> &msg, bool wait);

void delete_socket(zmq::socket_t *socket);


#endif
