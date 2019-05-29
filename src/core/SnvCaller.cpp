//
// Created by Chuanyi Zhang on 2019-05-23.
//

#include "SnvCaller.h"
#include <cmath>
#include <map>
#include <numeric>

using namespace moss;

inline uint8_t set_difference(uint8_t a, uint8_t b) {
    return a & ~b;
}

unsigned count_1bits[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

// TODO: limiting number of tumor samples?
SnvCaller::SnvCaller(int n_tumor_sample, double mu, double stepSize) : n_tumor_sample(n_tumor_sample), mu(mu),
                                                                       stepSize(stepSize) {
    gridSize = static_cast<int>(1 / stepSize + 1);
    logNoisePriorComplement = log(1 - mu);
    logPriorZComplement = logNoisePriorComplement - log((1 << n_tumor_sample) - 1);
    logMu = log1p(mu);
    logUniform = -log(gridSize - 1);
    eps = 0.1;
}

uint8_t SnvCaller::normal_calling(std::vector<Read> column, uint8_t ref) {
    std::map<uint8_t, int> count;
    for (const auto &r : column) {
        count[r.base]++;
    }
    for (const auto &c : count) {
        if ((c.first & ref) != 0) {
            if (abs(c.second / static_cast<double>(column.size()) - 0.5) <= eps) {
                return c.first | ref;
            } else if (1 - eps <= c.second / static_cast<double>(column.size())) {
                return c.first;
            }
        }
    }
    return ref;
}

Array3D
SnvCaller::likelihood(std::vector<std::vector<Read>> aligned, std::unordered_set<uint8_t> normal_bases,
                      std::unordered_set<uint8_t> tumor_base) {
    // pre-calculate
    auto n_gt = tumor_base.size();
    auto p_err = new double *[n_tumor_sample];
    auto is_normal = new bool *[n_tumor_sample];
    auto is_tumor = new bool *[n_tumor_sample];
    auto n_tumor = new int[n_tumor_sample * n_gt]{};
    auto n_normal = new int[n_tumor_sample]{};
    int idx_sample = 0;
    for (const auto &sample : aligned) {
        p_err[idx_sample] = new double[sample.size()];
        is_normal[idx_sample] = new bool[sample.size()];
        is_tumor[idx_sample] = new bool[sample.size() * n_gt];
        int idx_read = 0;
        for (auto r = sample.begin(); r != sample.end(); ++r, ++idx_read) {
            p_err[idx_sample][idx_read] = qphred2prob(r->qual);
            is_normal[idx_sample][idx_read] = normal_bases.count(r->base) > 0;
            n_normal[idx_sample] += is_normal[idx_sample][idx_read] ? 1 : 0;
            int idx_base = 0;
            for (const auto &tumorBase : tumor_base) {
                bool eq = r->base == tumorBase;
                is_tumor[idx_sample][idx_read * n_gt + idx_base] = eq;
                n_tumor[idx_sample * n_gt + idx_base] += eq ? 1 : 0;
                idx_base++;
            }
        }
        idx_sample++;
    }
    // n_tumor_sample x n_tumor_base x gridSize
    Array3D loglikelihood_3d;
    loglikelihood_3d.reserve(n_tumor_sample);
    idx_sample = 0;
    for (const auto &sample : aligned) {
        std::vector<std::vector<double>> sample_llh_2d;
        sample_llh_2d.reserve(tumor_base.size());
        int idx_base = 0;
        for (const auto &base : tumor_base) {
            std::vector<double> base_llh_1d;
            base_llh_1d.reserve(gridSize);
            for (int i = 0; i < gridSize; ++i) {
                double lhood = 0;
                double f = i * stepSize;
                unsigned long sample_size = sample.size();
                for (int j = 0; j < sample_size; ++j) {
                    if (is_normal[idx_sample][j]) {
                        lhood += log(f * p_err[idx_sample][j] / 3 + (1 - f) * (1 - p_err[idx_sample][j]));
                    } else if (is_tumor[idx_sample][j * n_gt + idx_base]) {
                        lhood += log(f * (1 - p_err[idx_sample][j]) + (1 - f) * p_err[idx_sample][j] / 3);
                    } else {
                        // f * err / 3 + (1-f) * err / 3
                        lhood += log(p_err[idx_sample][j]) - log(3);
                    }
                }
                unsigned int n_t = n_tumor[idx_sample * n_gt + idx_base];
                unsigned long n_n = n_normal[idx_sample];
                base_llh_1d.push_back(lhood + log(trinomial(sample_size - n_n - n_t, n_n, n_t)));
            }
            sample_llh_2d.emplace_back(base_llh_1d);
            ++idx_base;
        }
        loglikelihood_3d.emplace_back(sample_llh_2d);
        ++idx_sample;
    }

    // free
    for (int idx_sample = 0; idx_sample < n_tumor_sample; ++idx_sample) {
        delete[](p_err[idx_sample]);
        delete[](is_normal[idx_sample]);
        delete[](is_tumor[idx_sample]);
    }
    delete[](p_err);
    delete[](is_normal);
    delete[](is_tumor);
    delete[](n_tumor);
    delete[](n_normal);
    return loglikelihood_3d;
}

Array3D SnvCaller::likelihood(std::vector<std::vector<Read>> aligned, uint8_t normal_bases, uint8_t tumor_base) {
    std::unordered_set<uint8_t> normal, tumor;
    for (int i = 0; i < 4; ++i) {
        if (normal_bases & (1 << i)) {
            normal.insert(static_cast<uint8_t >(1 << i));
        }
        if (tumor_base & (1 << i)) {
            tumor.insert(static_cast<uint8_t >(1 << i));
        }
    }
    return SnvCaller::likelihood(aligned, normal, tumor);
}

double SnvCaller::calling(Pileups pile, uint8_t &normal_gt, uint8_t &tumor_gt) {
    const std::vector<std::vector<Read>> &columns = pile.get_read_columns();
    uint8_t ref = pile.get_ref();
    normal_gt = normal_calling(columns[0], ref);
    uint8_t tumor_bases = set_difference(0x0f, ref);
    Array3D lhood = likelihood(std::vector<std::vector<moss::Read>>(columns.begin() + 1, columns.end()), ref,
                               tumor_bases);
    unsigned int n_tumor_bases = count_1bits[tumor_bases];

    // pre-compute integral of log likelihood under z = 0, 1
    // integral over frequency, P(D_i | Z_i=z_i, Gn)
    std::vector<std::vector<double> > llh_integral(n_tumor_sample);
    int idx_sample = 0;
    for (const auto &sample : lhood) {
        int idx_base = 0;
        llh_integral[idx_sample].resize(n_tumor_bases);
        for (const auto &sample_base : sample) {
            llh_integral[idx_sample][idx_base] = log_sum_exp(sample_base) + logUniform;
            idx_base++;
        }
        idx_sample++;
    }
    // sum over 2^m-1 of z, and tumor nucleotide
    double nume = 0,
        deno = 0,
        max_tumor_evidence = 0;
    int tumor_gt_idx;
    for (int idx_nuc = 0; idx_nuc < n_tumor_bases; ++idx_nuc) {
        double evidence = 0;
        for (int z = 0; z < (1 << n_tumor_sample); ++z) {
            double llh = 0;
            for (int idx_sample = 0; idx_sample < n_tumor_sample; ++idx_sample) {
                int z_sample = (z >> idx_sample) & 1;
                llh += z_sample ? llh_integral[idx_sample][idx_nuc] : lhood[idx_sample][idx_nuc][0];
            }
            if (z == 0) {
                double temp = exp(llh + logMu);
                nume += temp;
                evidence += temp;
            } else {
                evidence += exp(llh + logPriorZComplement);
            }
            if (max_tumor_evidence <= evidence) {
                max_tumor_evidence = evidence;
                // TODO: tumor_gt =
                tumor_gt_idx = idx_nuc;
            }
            deno += evidence;
        }
    }
    double log_prob_non_soma;
    if (nume > 0) {
        log_prob_non_soma = log(nume) - log(deno);
    } else {
        log_prob_non_soma = float('-inf');
    }
    return log_prob_non_soma;
}

double moss::qphred2prob(int qphred) {
    return pow(10.f, -static_cast<double>(qphred) / 10.f);
}

double moss::binom(unsigned int n, unsigned int k) {
    assert(n >= k);
    if (n == 0 || n == 1 || k == 0 || k == n) {
        return 1;
    }
    if (k == 1) {
        return n;
    }
    if (2 * k > n) {
        k = n - k;
    }
    double bin = n - k + 1;
    for (int i = 2, j = n - k + 2; i <= k; ++i, ++j) {
        bin *= j;
        bin /= i;
    }
    return bin;
}

double moss::trinomial(unsigned int s, unsigned int k, unsigned int t) {
    return binom(s + k + t, t) * binom(s + k, k);
}

template<typename T>
T moss::log_sum_exp(std::vector<T> array) {
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
