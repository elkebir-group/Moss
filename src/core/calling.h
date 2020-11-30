//
// Created by Chuanyi Zhang on 2019-05-23.
//

#ifndef MOSS_CALLING_H
#define MOSS_CALLING_H

#include <vector>
#include <limits>
#include <cmath>
#include "types.h"
#include "../io/vcf_io.h"

namespace moss {
    using Array3D = std::vector<std::vector<std::vector<double> > >;

    double qphred2prob(int qphred);

    /**
     * Log trinomial coefficient.
     * @return (s+k+t)! / s! / k! / t!.
     */
    double log_trinomial(unsigned long s, unsigned long k, unsigned long t);

    /**
     * perform log sum exp trick in one pass
     * @tparam T : numerical that supports add and log
     * @param array : array of log values
     * @return log summation
     */
    template<typename T>
    T log_sum_exp(std::vector<T> array) {
        T max_elem = -std::numeric_limits<T>::infinity(),
            accum{};
        for (const auto &item : array) {
            if (item >= max_elem) {
                accum *= exp(max_elem - item);
                accum += 1.f;
                max_elem = item;
            } else {
                accum += exp(item - max_elem);
            }
        }
        return max_elem + log(accum);
    }

    template<typename T>
    T log_sum_exp(T a, T b) {
        if (a >= b) {
            return a + log(exp(b - a) + 1.0);
        } else {
            return b + log(exp(a - b) + 1.0);
        }
    }

    template<typename T>
    T log_sum_exp_array(std::vector<T> array1, std::vector<T> array2) {
        T max_elem = -std::numeric_limits<T>::infinity(),
            accum{};
        for (int idx = 0; idx < array1.size(); idx++) {
            T item = array1[idx] + array2[idx];
            if (item >= max_elem) {
                accum *= exp(max_elem - item);
                accum += 1.f;
                max_elem = item;
            } else {
                accum += exp(item - max_elem);
            }
        }
        return max_elem + log(accum);
    }

    /**
     * Initialize max_elem to be -infinity.
     * @tparam T is the type of data.
     * @param max_elem is the maximum element so far.
     * @param accum is the accumulated log sum.
     */
    template<typename T>
    void log_sum_exp_init(T &max_elem, T &accum) {
        max_elem = -std::numeric_limits<T>::infinity();
        accum = T{};
    }

    /**
     * Single iteration of log sum exp trick.
     * @tparam T is the type of data.
     * @param max_elem is the maximum element so far.
     * @param accum is the accumulated log sum.
     * @param item is the new element.
     */
    template<typename T>
    void log_sum_exp_iter(T &max_elem, T &accum, T item) {
        if (item >= max_elem) {
            accum *= exp(max_elem - item);
            accum += 1.f;
            max_elem = item;
        } else {
            accum += exp(item - max_elem);
        }
    }

    /**
     * Finalize the log sum exp trick.
     * @tparam T is the type of data.
     * @param max_elem is the maximum element so far.
     * @param accum is the accumulated log sum.
     * @return The total log(sum(exp())).
     */
    template<typename T>
    T log_sum_exp_final(T &max_elem, T &accum) {
        return max_elem + log(accum);
    }

    class SnvCaller {
    private:
        VcfReader<RecData> normal_result;
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
        bool is_ignore0;

        std::vector<std::vector<double>> p_err;
        std::vector<std::vector<char>> is_normal;
        std::vector<std::vector<char>> is_tumor;
        std::vector<bool> is_empty;
        std::vector<bool> is_all_ref;
        std::vector<int> n_tumor;   //* count of tumor allels in samples: n_tumor_sample x n_tumor_bases
        std::vector<int> n_normal;  //* count of normal allels in samples: n_tumor_sample
        std::vector<double> p_err_normal;
        std::vector<int> is_normal_normal;
        std::vector<char> is_tumor_normal;
        std::vector<int> n_tumor_in_normal;
        int n_normal_in_normal;
        std::vector<std::vector<double>> likelihoods_normal;
        std::vector<double> log_TIN;
        std::vector<double> log_beta_pdf;

        BaseSet normal_calling(const std::vector<Read> &column, uint8_t ref) const;

        BaseSet normal_calling(const std::string &contig, locus_t pos, uint8_t ref);

        Array3D likelihoods;        //* n_tumor_sample x n_tumor_base x gridSize

        /**
         * Calculate log likelihood of each samples given tumor base and VAF,
         * P( D_j | f_j, z_j, g_t, g_n).
         * modify member lhood, a 3D vector of size n_tumor_sample x n_tumor_base x gridSize
         * @param aligned are aligned columns of tumor samples.
         * @param normal_bases is a set of normal bases.
         * @param tumor_base is the available tumor bases.
         */
        void calc_likelihood(const std::vector<std::vector<Read>> &aligned, BaseSet normal_bases, BaseSet tumor_base);

    public:
        /**
         * Construct a SNV caller
         *
         * @param n_tumor_sample is the number of tumor samples.
         * @param normal is the name of the normal VCF file.
         * @param is_ignore0 indicates whether ignore samples with all reads same as normal allele.
         * @param mu is the parameter for mutational rate.
         * @param max_depth is the advised max-depth in data, will automatically increase when reached.
         * @param grid_size is the number of grid points approximating the integration.
         */
        SnvCaller(int n_tumor_sample, const std::string& normal, bool is_ignore0, double mu = 1.0 - 5e-4,
                   int max_depth = 500, int grid_size = 21);

        ~SnvCaller();

        /**
         * Calculate the somatic SNV probability and fill the annotations.
         * @param chrom is the chromosome.
         * @param pos is the position.
         * @param pile is the column in pileup.
         * @param anno is the annotation to be filled.
         * @return The somatic SNV probability in log space.
         */
        double
        calling(const std::string &chrom, locus_t pos, const Pileups &pile, Annotation &anno);
    };
}


#endif //MOSS_CALLING_H
