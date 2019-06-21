#include "utils/string_util.h"


std::string &ltrim(std::string &s) {
    auto it = s.begin();
    for (; it != s.end(); ++it) {
        if (!std::isspace(*it)) {
            break;
        }
    }
    s.erase(s.begin(), it);
    return s;
}

// trim from end
std::string &rtrim(std::string &s) {
    auto it = s.end() - 1;
    for (; it != s.begin() - 1; --it) {
        if (!std::isspace(*it)) {
            break;
        }
    }
    s.erase(it + 1, s.end());
    return s;
}

std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}
