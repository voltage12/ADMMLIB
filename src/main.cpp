#include <iostream>
#include <string>
#include <unistd.h>

#include "utils/properties.h"
#include "coordinator.h"


int main(int argc, char **argv) {
    if (argc < 3) {
        std::cout << "program hostfile_path configuration_file_path" << std::endl;
    } else {
        Properties properties(argv[2]);
        Coordinator coordinator(argv[1], properties);
        coordinator.run();
    }
    return 0;
}