#include <iostream>
#include <fstream>
#include <sstream>

#include "optimizer/sparse_sample_set.h"


SparseSampleSet::SparseSampleSet(const std::string &data_file_path) {
    all_feature_num = 0;
    sample_num = 0;
    dim = 0;
    std::ifstream file_reader(data_file_path);
    if (file_reader.fail()) {
        std::cerr << "打开数据文件 " << data_file_path << " 失败" << std::endl;
        exit(-1);
    }
    int label;
    std::string temp_str;
    std::string line;
    std::stringstream line_reader;
    while (std::getline(file_reader, line)) {
        line.erase(0, line.find_first_not_of(' '));
        line.erase(line.find_last_not_of(' ') + 1);
        if (line.empty()) {
            std::cerr << "出现空行" << std::endl;
            continue;
        } else {
            line_reader.clear();
            line_reader.str(line);
            if (!(line_reader >> label) || !(label == 1 || label == -1 || label == 0)) {
                std::cerr << "读取标签失败" << std::endl;
                exit(-1);
            }
            while (line_reader >> temp_str) {
                size_t pos = temp_str.find_first_of(':');
                if (pos == std::string::npos) {
                    std::cerr << "读取数据文件 " << data_file_path << " 时出错，第" << sample_num + 1 << "行缺少:符号" << std::endl;
                    exit(-1);
                } else if (pos == 0) {
                    std::cerr << "读取数据文件 " << data_file_path << " 时出错，第" << sample_num + 1 << "行缺少index" << std::endl;
                    exit(-1);
                } else if (pos == temp_str.size() - 1) {
                    std::cerr << "读取数据文件 " << data_file_path << " 时出错，第" << sample_num + 1 << "行缺少value" << std::endl;
                    exit(-1);
                }
                std::string index = temp_str.substr(0, pos);
                std::string value = temp_str.substr(pos + 1);
                size_t idx;
                std::stoi(index, &idx);
                if (idx != index.size()) {
                    std::cerr << "读取数据文件 " << data_file_path << " 时出错，第" << sample_num + 1 << "行index非法" << std::endl;
                    exit(-1);
                }
                std::stod(value, &idx);
                if (idx != value.size()) {
                    std::cerr << "读取数据文件 " << data_file_path << " 时出错，第" << sample_num + 1 << "行value非法" << std::endl;
                    exit(-1);
                }
                ++all_feature_num;
            }
            //每条样本以一个index=-1的Feature表示结束
            ++all_feature_num;
            ++sample_num;
        }
    }

    label_list = new int[sample_num];
    sample_list = new Feature *[sample_num];
    sample_space = new Feature[all_feature_num];

    //让文件读取指针归位
    file_reader.clear();
    file_reader.seekg(0, std::ifstream::beg);
    int i = 0, j = 0;
    while (std::getline(file_reader, line)) {
        line.erase(0, line.find_first_not_of(' '));
        line.erase(line.find_last_not_of(' ') + 1);
        if (line.empty()) {
            continue;
        } else {
            //一条样本的起始Feature的地址
            sample_list[i] = &sample_space[j];
            line_reader.clear();
            line_reader.str(line);
            line_reader >> label_list[i];
            if (label_list[i] != 1) {
                label_list[i] = -1;
            }
            while (line_reader >> temp_str) {
                size_t pos = temp_str.find_first_of(':');
                std::string index = temp_str.substr(0, pos);
                std::string value = temp_str.substr(pos + 1);
                sample_space[j].index = atoi(index.c_str());
                sample_space[j].value = atof(value.c_str());
                if (sample_space[j].index > dim) {
                    dim = sample_space[j].index;
                }
                ++j;
            }
            //一条样本以index=-1的Feature结尾
            sample_space[j++].index = -1;
            ++i;
        }
    }
    if (i != sample_num || j != all_feature_num) {
        std::cerr << "读取数据出错" << std::endl;
        exit(-1);
    }
    file_reader.close();
}

SparseSampleSet::~SparseSampleSet() {
    delete[] label_list;
    delete[] sample_list;
    delete[] sample_space;
}

int SparseSampleSet::get_label(int n) {
    if (n < sample_num) {
        return label_list[n];
    } else {
        std::cerr << "越界，最大样本数量为" << sample_num << std::endl;
        exit(-1);
    }
}

const Feature *SparseSampleSet::get_sample(int n) {
    if (n < sample_num) {
        return sample_list[n];
    } else {
        std::cerr << "越界，最大样本数量为" << sample_num << std::endl;
        exit(-1);
    }
}

double SparseSampleSet::dot(int n, const double *x) {
    double ret = 0;
    const Feature *sample = get_sample(n);
    while (sample->index != -1) {
        ret += x[(sample->index - 1)] * sample->value;
        ++sample;
    }
    return ret;
}