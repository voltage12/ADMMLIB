#include <iostream>
#include <sstream>
#include <fstream>

#include "utils/properties.h"
#include "utils/string_util.h"


Properties::Properties(const std::string &filename) {
    std::ifstream reader(filename);
    if (reader.fail()) {
        std::cerr << "无法打开配置文件:" << filename << std::endl;
        exit(-1);
    }

    std::string line;
    while (std::getline(reader, line)) {
        std::size_t pos = line.find_first_of('#');
        if (pos != std::string::npos) {
            line.erase(pos);
        }
        trim(line);
        if (line.empty()) {
            continue;
        }
        pos = line.find_first_of(':');
        if (pos == std::string::npos) {
            std::cerr << "配置名和值之间应该用':'分开 " << line << std::endl;
            exit(-1);
        }
        std::string key = line.substr(0, pos);
        std::string value = line.substr(pos + 1);
        properties[key] = value;
    }
    reader.close();
}

std::string Properties::get_string(const std::string &property_name) {
    return properties.at(property_name);
}

int Properties::get_int(const std::string &property_name) {
    std::stringstream type_converter(properties.at(property_name));
    int value = 0;
    if (!(type_converter >> value)) {
        std::cerr << property_name << "应该为一个整数" << std::endl;
        exit(-1);
    }
    return value;
}

double Properties::get_double(const std::string &property_name) {
    std::stringstream type_converter(properties.at(property_name));
    double value = 0;
    if (!(type_converter >> value)) {
        std::cerr << property_name << "应该为一个浮点数" << std::endl;
        exit(-1);
    }
    return value;
}

bool Properties::get_bool(const std::string &property_name) {
    if (properties.at(property_name) == "true") {
        return true;
    } else if (properties.at(property_name) == "false") {
        return false;
    } else {
        std::cerr << property_name << "的值必须为true或者false" << std::endl;
        exit(-1);
    }
}

bool Properties::has_property(const std::string &property_name) {
    return properties.count(property_name) != 0;
}

void Properties::check_property(const std::string &property_name) {
    if (!has_property(property_name)) {
        std::cerr << "配置文件中并没有参数" << property_name << std::endl;
        exit(-1);
    }
}

void Properties::print() {
    std::cout << "#**************Configuration**************" << std::endl;
    for (auto it = properties.begin(); it != properties.end(); ++it) {
        std::cout << it->first << ": " << it->second << std::endl;
    }
    std::cout << "#*****************************************" << std::endl;
}


