#ifndef PROPERTIES_H
#define PROPERTIES_H

#include <map>
#include <string>


class Properties {
public:
    explicit Properties(const std::string &filename);

    std::string get_string(const std::string &property_name);

    int get_int(const std::string &property_name);

    double get_double(const std::string &property_name);

    bool get_bool(const std::string &property_name);

    bool has_property(const std::string &property_name);

    void check_property(const std::string &property_name);

    void print();

private:
    std::map<std::string, std::string> properties;
};


#endif
