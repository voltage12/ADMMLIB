#ifndef SPARSE_SAMPLE_SET_H
#define SPARSE_SAMPLE_SET_H

#include <string>


struct Feature {
    int index;
    double value;
};

class SparseSampleSet {
public:
    SparseSampleSet(const std::string &data_file_path);

    ~SparseSampleSet();

    inline int get_dim() { return dim; }

    inline int get_sample_num() { return sample_num; }

    const Feature *get_sample(int n);

    double dot(int n, const double *x);

    int get_label(int n);

private:
    int sample_num; /* 样本数量 */
    int dim; /* 样本维度 */
    int all_feature_num;
    int *label_list;

    Feature **sample_list;
    Feature *sample_space;
};


#endif
