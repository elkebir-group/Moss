//
// Created by Chuanyi Zhang on 2019-05-23.
//

#ifndef MOSS_CALLING_H
#define MOSS_CALLING_H

#include <vector>
#include "types.h"
#include "../io/vcf_io.h"

namespace moss {
    using Array3D = std::vector<std::vector<std::vector<double> > >;

    double qphred2prob(int qphred);

    double log_trinomial(unsigned long s, unsigned long k, unsigned long t);

    /*!
     * perform log sum exp trick in one pass
     * @tparam T : numerical that supports add and log
     * @param array : array of log values
     * @return log summation
     */
    template<typename T>
    T log_sum_exp(std::vector<T> array);

    template<typename T>
    T log_sum_exp(T a, T b);

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
        VcfReader normal_result;
        int n_tumor_sample,
            gridSize,
            max_depth;
        // neucleotide = {'T', 'A', 'G', 'C'};
        // genotype = {"ref": 0, "het": 0.5, "hom": 1};
        double stepSize,
            logUniform,
            mu,
            *logPriorZComplement,
            *logAll0,
            logMu,
            logNoisePriorComplement,
            eps;

        std::vector<std::vector<double>> p_err;
        std::vector<std::vector<char>> is_normal;
        std::vector<std::vector<char>> is_tumor;
        std::vector<bool> is_empty;
        std::vector<int> n_tumor;   //* count of tumor allels in samples: n_tumor_sample x n_tumor_bases
        std::vector<int> n_normal;  //* count of normal allels in samples: n_tumor_sample
        std::vector<double> p_err_normal;
        std::vector<int> is_normal_normal;
        std::vector<int> is_tumor_normal;
        std::vector<double> likelihoods_normal;

        BaseSet normal_calling(const std::vector<Read> &column, uint8_t ref);

        BaseSet normal_calling(const std::string &contig, locus_t pos, uint8_t ref);

        Array3D likelihoods;        //* n_tumor_sample x n_tumor_base x gridSize

        /*!
         * Calculate log likelihood of each samples given tumor base and VAF,
         * P( D_j | f_j, z_j, g_t, g_n).
         * modify member lhood, a 3D vector of size n_tumor_sample x n_tumor_base x gridSize
         * @param aligned : aligned columns of tumor samples
         * @param normal_bases : base set of normal given by normal_calling
         * @param tumor_base : return MLE tumor base
         */
        void calc_likelihood(const std::vector<std::vector<Read>> &aligned, BaseSet normal_bases, BaseSet tumor_base);
        double in_normal(const Pileups &pile, BaseSet &normal_gt, const uint8_t &tumor_gt);

    public:

        SnvCaller(int n_tumor_sample, const std::string& normal, double mu = 1.0 - 5e-4, int max_depth = 500,
                  int grid_size = 21);

        ~SnvCaller();

        double
        calling(const std::string &chrom, locus_t pos, const Pileups &pile, unsigned long &Z, Annotation &anno);
    };
}


#endif //MOSS_CALLING_H
