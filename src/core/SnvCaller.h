//
// Created by Chuanyi Zhang on 2019-05-23.
//

#ifndef MOSS_SNVCALLER_H
#define MOSS_SNVCALLER_H

#include <vector>
#include <unordered_set>
#include "types.h"

namespace moss {
    using Array3D = std::vector<std::vector<std::vector<double> > >;

    double qphred2prob(int qphred);

    double binom(unsigned int n, unsigned int k);

    double trinomial(unsigned int s, unsigned int k, unsigned int t);

    template<typename T>
    T log_sum_exp(std::vector<T> array);

    class SnvCaller {
    private:
        int n_tumor_sample,
                gridSize;
        // neucleotide = {'T', 'A', 'G', 'C'};
        // genotype = {"ref": 0, "het": 0.5, "hom": 1};
        double stepSize,
                logUniform,
                mu,
                logPriorZComplement,
                logMu,
                logNoisePriorComplement,
                eps;

        uint8_t normal_calling(std::vector<Read> column, uint8_t ref);

        Array3D likelihood(std::vector<std::vector<Read>> aligned, std::unordered_set<uint8_t> normal_bases,
                           std::unordered_set<uint8_t> tumor_base);

    public:

        SnvCaller(int n_tumor_sample, double mu = 1.0 - 5e-6, double stepSize = 0.05);

        Array3D likelihood(std::vector<std::vector<Read>> aligned, uint8_t normal_bases, uint8_t tumor_base);

        double calling(Pileups pile, uint8_t &normal_gt, uint8_t &tumor_gt);
    };
}


#endif //MOSS_SNVCALLER_H
