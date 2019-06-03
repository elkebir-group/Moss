//
// Created by Chuanyi Zhang on 2019-05-23.
//

#ifndef MOSS_SNVCALLER_H
#define MOSS_SNVCALLER_H

#include <vector>
#include "types.h"

namespace moss {
    using Array3D = std::vector<std::vector<std::vector<double> > >;

    double qphred2prob(int qphred);

    double binom(unsigned int n, unsigned int k);

    double trinomial(unsigned long s, unsigned long k, unsigned long t);

    /*!
     * perform log sum exp trick in one pass
     * @tparam T : numerical that supports add and log
     * @param array : array of log values
     * @return log summation
     */
    template<typename T>
    T log_sum_exp(std::vector<T> array);

    /*!
     * Initialize max_elem to be -infinity, return accumulate.
     * @tparam T
     * @param max_elem
     * @return accumulate
     */
    template<typename T>
    T &log_sum_exp_init(T &max_elem);

    /*!
     * Single iteration of log sum exp trick.
     * @tparam T
     * @param max_elem
     * @param accum
     * @param item
     */
    template<typename T>
    void log_sum_exp_iter(T &max_elem, T &accum, T item);

    template<typename T>
    T log_sum_exp_final(T &max_elem, T &accum);

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

        BaseSet normal_calling(std::vector<Read> column, uint8_t ref);

//        Array3D likelihood(std::vector<std::vector<Read>> aligned, std::unordered_set<uint8_t> normal_bases,
//                           std::unordered_set<uint8_t> tumor_base);

    public:

        SnvCaller(int n_tumor_sample, double mu = 1.0 - 5e-6, double stepSize = 0.05);

        Array3D likelihood(std::vector<std::vector<Read>> aligned, BaseSet normal_bases, BaseSet tumor_base);

        double calling(Pileups pile, BaseSet &normal_gt, uint8_t &tumor_gt);
    };
}


#endif //MOSS_SNVCALLER_H
