//
// Created by Chuanyi Zhang on 2019-05-23.
//

#ifndef MOSS_CALLING_H
#define MOSS_CALLING_H

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
     * @param accum
     * @return accumulate
     */
    template<typename T>
    void log_sum_exp_init(T &max_elem, T &accum);

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
            gridSize,
            max_depth;
        // neucleotide = {'T', 'A', 'G', 'C'};
        // genotype = {"ref": 0, "het": 0.5, "hom": 1};
        double stepSize,
                logUniform,
                mu,
                logPriorZComplement,
                logMu,
                logNoisePriorComplement,
                eps;

        double** p_err;
        bool** is_normal;
        bool** is_tumor;
        int* n_tumor;
        int* n_normal;

        BaseSet normal_calling(const std::vector<Read> &column, uint8_t ref);

        /*!
         * Calculate log likelihood of each samples given tumor base and VAF,
         * P( D_j | f_j, z_j, g_t, g_n).
         * Return a 3D vector of size n_tumor_sample x n_tumor_base x gridSize
         * @param aligned : aligned columns of tumor samples
         * @param normal_bases : base set of normal given by normal_calling
         * @param tumor_base : return MLE tumor base
         * @return : 3D vector
         */
        Array3D likelihood(const std::vector<std::vector<Read>> &aligned, BaseSet normal_bases, BaseSet tumor_base);

    public:

        SnvCaller(int n_tumor_sample, double mu = 1.0 - 5e-6, double stepSize = 0.05, int max_depth = 500);

        virtual ~SnvCaller();

        double calling(const Pileups &pile, BaseSet &normal_gt, uint8_t &tumor_gt, unsigned long &Z);
    };
}


#endif //MOSS_CALLING_H
