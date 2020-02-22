//
// Created by Chuanyi Zhang on 2019-05-23.
//

#include "calling.h"
#include <cmath>
#include <map>
#include <cassert>
#include <cstring>
#include <limits>

using namespace moss;

// TODO: limited number of tumor samples?
// - Use custom arbitrary long binary indicator?
SnvCaller::SnvCaller(int n_tumor_sample, const std::string& normal, bool is_ignore0, double mu, int max_depth,
                     int grid_size)
    : n_tumor_sample(n_tumor_sample),
      normal_result(VcfReader(normal)),
      mu(mu),
      max_depth(max_depth),
      gridSize(grid_size),
      is_ignore0(is_ignore0),
      is_empty(n_tumor_sample),
      is_all_ref(n_tumor_sample),
      n_tumor(n_tumor_sample * 3),
      n_normal(n_tumor_sample),
      p_err(n_tumor_sample, std::vector<double>(max_depth)),
      is_normal(n_tumor_sample, std::vector<char>(max_depth)),
      is_tumor(n_tumor_sample, std::vector<char>(max_depth * 3)),
      p_err_normal(max_depth),
      is_normal_normal(max_depth),
      is_tumor_normal(max_depth),
      likelihoods_normal(grid_size),
      likelihoods(n_tumor_sample, std::vector<std::vector<double>>(3, std::vector<double>(grid_size))) {
    stepSize = 1.0 / (gridSize - 1);
    logMu = log1p(mu - 1);
    logNoisePriorComplement = log(1 - mu);
    logPriorZComplement = new double[n_tumor_sample];
    logAll0 = new double[n_tumor_sample];
    for (int idx = 0; idx < n_tumor_sample; idx++) {
        logPriorZComplement[idx] = logNoisePriorComplement - log((1ULL << (idx + 1)) - 1);
        // log(1 - \mu * 2^m / (2^m - 1)) = log(2^m - 1 - 2^m(1-MU)) - log(2^m - 1)
        if (idx == 0) {
            logAll0[idx] = log1p(mu * (1ULL << (idx + 1)) - 2) - log1p((1ULL << (idx + 1)) - 2);
        } else {
            logAll0[idx] = log(mu * (1ULL << (idx + 1)) - 1) - log((1ULL << (idx + 1)) - 1);
        }
    }
    logUniform = -log(gridSize);
    eps = 0.1;
}

SnvCaller::~SnvCaller() {
    delete[](logPriorZComplement);
}

BaseSet SnvCaller::normal_calling(const std::vector<Read> &column, uint8_t ref) {
    std::map<uint8_t, int> count;
    for (const auto &r : column) {
        count[r.base]++;
    }
    for (const auto &c : count) {
        if ((c.first & ref) != 0) {
            if (abs(c.second / static_cast<double>(column.size()) - 0.5) <= eps) {
                return BaseSet(c.first | ref);
            } else if (1 - eps <= c.second / static_cast<double>(column.size())) {
                return BaseSet(c.first);
            }
        }
    }
    return BaseSet(ref);
}

BaseSet SnvCaller::normal_calling(const std::string &contig, locus_t pos, uint8_t ref) {
    RecData search = normal_result.find(contig, pos);
    if (search.bases.is_valid() && search.is_pass) {
        return search.bases;
    } else {
        return BaseSet(ref);
    }
}

void
SnvCaller::calc_likelihood(const std::vector<std::vector<Read>> &aligned, BaseSet normal_bases, BaseSet tumor_base) {
    // pre-calculate
    auto n_gt = tumor_base.size();
    std::fill(n_normal.begin(), n_normal.end(), 0);
    std::fill(n_tumor.begin(), n_tumor.end(), 0);
    int idx_sample = 0;
    for (const auto &sample : aligned) {
        int idx_read = 0;
        for (auto r = sample.begin(); r != sample.end(); ++r, ++idx_read) {
            if (idx_read >= max_depth) {
                max_depth *= 2;
                for (int i = 0; i < n_tumor_sample; ++i) {
                    p_err[i].resize(max_depth);
                    is_normal[i].resize(max_depth);
                    is_tumor[i].resize(max_depth * 3);
                }
                p_err_normal.resize(max_depth);
                is_normal_normal.resize(max_depth);
                is_tumor_normal.resize(max_depth);
            }
            p_err[idx_sample].at(idx_read) = qphred2prob(r->qual);
            is_normal[idx_sample].at(idx_read) = normal_bases.contain(r->base);
            n_normal[idx_sample] += is_normal[idx_sample][idx_read] ? 1 : 0;
            int idx_base = 0;
            for (const auto &tumorBase : tumor_base) {
                bool eq = r->base == tumorBase;
                is_tumor[idx_sample].at(idx_read * n_gt + idx_base) = eq;
                n_tumor[idx_sample * n_gt + idx_base] += eq ? 1 : 0;
                idx_base++;
            }
        }
        idx_sample++;
    }
    // n_tumor_sample x n_tumor_base x gridSize
    idx_sample = 0;
    for (const auto &sample : aligned) {
        int idx_base = 0;
        for (const auto &base : tumor_base) {
            for (int idx_step = 0; idx_step < gridSize; ++idx_step) {
                double lhood = 0;
                double f = idx_step * stepSize;
                unsigned long sample_size = sample.size();
                assert(sample_size < max_depth);
                for (int j = 0; j < sample_size; ++j) {
                    if (is_normal[idx_sample][j]) {
                        lhood += log(f * p_err[idx_sample][j] / 3 +
                                     (1 - f) * (1 + p_err[idx_sample][j] * (double(normal_bases.size()) - 4.0) / 3));
                        assert(lhood <= 0);
                    } else if (is_tumor[idx_sample][j * n_gt + idx_base]) {
                        lhood += log(f * (1 - p_err[idx_sample][j]) + (1 - f) * p_err[idx_sample][j] / 3);
                        assert(lhood <= 0);
                    } else {
                        // f * err / 3 + (1-f) * err / 3
                        lhood += log(p_err[idx_sample][j]) - log(3);
                        assert(lhood <= 0);
                    }
                }
                unsigned int n_t = n_tumor[idx_sample * n_gt + idx_base];
                unsigned long n_n = n_normal[idx_sample];
                double coeff = log_trinomial(sample_size - n_n - n_t, n_n, n_t);
                assert(!std::isinf(coeff));
                likelihoods[idx_sample][idx_base][idx_step] = lhood + coeff;
            }
            ++idx_base;
        }
        ++idx_sample;
    }
}

double SnvCaller::in_normal(const Pileups &pile, BaseSet &normal_gt, const uint8_t &tumor_gt) {
    const std::vector<std::vector<Read>> &columns = pile.get_read_columns();
    auto aligned = columns[0];
    double log_in_normal{};
    // pre-calculate
    int n_normal = 0,
        n_tumor = 0;
    int idx_read = 0;
    for (auto r = aligned.begin(); r != aligned.end(); ++r, ++idx_read) {
        if (idx_read >= max_depth) {
            max_depth *= 2;
            for (int i = 0; i < n_tumor_sample; ++i) {
                p_err[i].resize(max_depth);
                is_normal[i].resize(max_depth);
                is_tumor[i].resize(max_depth * 3);
            }
            p_err_normal.resize(max_depth);
            is_normal_normal.resize(max_depth);
            is_tumor_normal.resize(max_depth);
        }
        p_err_normal.at(idx_read) = qphred2prob(r->qual);
        is_normal_normal.at(idx_read) = normal_gt.contain(r->base);
        n_normal += is_normal_normal[idx_read] ? 1 : 0;

        bool eq = r->base == tumor_gt;
        is_tumor_normal.at(idx_read) = eq;
        n_tumor += eq ? 1 : 0;
    }
    for (int idx_step = 0; idx_step < gridSize; ++idx_step) {
        double lhood = 0;
        double f = idx_step * stepSize;
        unsigned long sample_size = aligned.size();
        assert(sample_size < max_depth);
        for (int j = 0; j < sample_size; ++j) {
            if (is_normal_normal[j]) {
                lhood += log(f * p_err_normal[j] / 3 +
                                (1 - f) * (1 + p_err_normal[j] * (double(normal_gt.size()) - 4.0) / 3));
                assert(lhood <= 0);
            } else if (is_tumor_normal[j]) {
                lhood += log(f * (1 - p_err_normal[j]) + (1 - f) * p_err_normal[j] / 3);
                assert(lhood <= 0);
            } else {
                // f * err / 3 + (1-f) * err / 3
                lhood += log(p_err_normal[j]) - log(3);
                assert(lhood <= 0);
            }
        }
        double coeff = log_trinomial(sample_size - n_normal - n_tumor, n_normal, n_tumor);
        assert(!std::isinf(coeff));
        likelihoods_normal[idx_step] = lhood + coeff;
    }
    double llh_integral{};
    llh_integral = log_sum_exp(likelihoods_normal) + logUniform;
    log_in_normal = (likelihoods_normal[0] + logMu) - (log_sum_exp(llh_integral + logNoisePriorComplement, likelihoods_normal[0] + logMu));
    return log_in_normal;
}

double
SnvCaller::calling(const std::string &chrom, locus_t pos, const Pileups &pile, unsigned long &Z, Annotation &annos) {
    annos.tumor_gt = uint8_t(IUPAC_nuc::EQ);
    const std::vector<std::vector<Read>> &columns = pile.get_read_columns();
    uint8_t ref = pile.get_ref();
    if (normal_result.empty()) {
        annos.normal_gt = normal_calling(columns[0], ref);
    } else {
        annos.normal_gt = normal_calling(chrom, pos, ref);
    }
    BaseSet tumor_bases = annos.normal_gt.complement();
    calc_likelihood(std::vector<std::vector<moss::Read>>(columns.begin() + 1, columns.end()), annos.normal_gt, tumor_bases);
    size_t n_potential_tumor_sample = 0;     // n_potential_tumor_sample means: is_ignore0 ? not all=ref : not empty
    annos.cnt_read[0] = columns[0].size();
    for (unsigned long i = 1; i < columns.size(); i++) {
        if (columns[i].empty()) {
            is_empty[i - 1] = true;
        } else {
            is_empty[i - 1] = false;
            if (is_ignore0) {
                is_all_ref[i-1] = true;
                for (auto &&base : columns[i]) {
                    if (base.base != ref) {
                        is_all_ref[i-1] = false;
                        n_potential_tumor_sample++;
                        break;
                    }
                }
            } else {
                n_potential_tumor_sample++;
            }
        }
        annos.cnt_read[i] = columns[i].size();
    }
    if (n_potential_tumor_sample == 0) {
        Z = 0;
        annos.tumor_gt = uint8_t(IUPAC_nuc::EQ);
        annos.quality = -0;
        annos.log_t_in_normal = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n_tumor_sample + 1; i++) {
            annos.cnt_tumor[i] = 0;
            annos.genotype[i] = 0;
            annos.zq[i] = 0;
        }
        std::fill(annos.cnt_type_strand.begin(), annos.cnt_type_strand.end(), 0);
        return 0;
    }

    // pre-compute integral of log likelihood under z = 0, 1
    // integral over frequency, P(D_i | Z_i=z_i, Gn)
    std::vector<std::vector<double> > llh_integral(n_tumor_sample);
    {int idx_sample = 0;
    for (const auto &sample : likelihoods) {
        llh_integral[idx_sample].resize(tumor_bases.size());
        for (int idx_base = 0; idx_base < tumor_bases.size(); ++idx_base) {
            const auto &sample_base = sample[idx_base];
            llh_integral[idx_sample][idx_base] = log_sum_exp(sample_base) + logUniform;
        }
        idx_sample++;
    }}
    // sum over 2^m-1 of z, and tumor nucleotide
    double max_nume_elem,
        nume,
        max_deno_elem,
        deno,
        max_tumor_evidence = -std::numeric_limits<double>::infinity();

    log_sum_exp_init(max_nume_elem, nume);
    log_sum_exp_init(max_deno_elem, deno);
    int tumor_gt_idx;
    int idx_nuc = 0;
    for (const auto &tumor_base : tumor_bases) {
        double nume_gt{},
            evidence{};
        for (int idx_sample = 0; idx_sample < n_tumor_sample; ++idx_sample) {
            bool skip;
            if (is_ignore0)
                skip = is_empty[idx_sample] || is_all_ref[idx_sample];
            else
                skip = is_empty[idx_sample];
            if (skip) {
                continue;
            }
            nume_gt += likelihoods[idx_sample][idx_nuc][0];
            evidence += log_sum_exp(llh_integral[idx_sample][idx_nuc], likelihoods[idx_sample][idx_nuc][0]);
        }
        evidence = log_sum_exp(evidence + logPriorZComplement[n_potential_tumor_sample - 1],
                               nume_gt + logAll0[n_potential_tumor_sample - 1]);
        if (max_tumor_evidence <= evidence) {
            max_tumor_evidence = evidence;
            annos.tumor_gt = tumor_base;
            tumor_gt_idx = idx_nuc;
        }
        log_sum_exp_iter(max_nume_elem, nume, nume_gt + logMu);
        log_sum_exp_iter(max_deno_elem, deno, evidence);
        idx_nuc++;
    }

    annos.log_t_in_normal = -10 * in_normal(pile, annos.normal_gt, annos.tumor_gt);

    annos.genotype[0] = 0;
    int normal_tumor_allele_count = 0;
    for (auto &item : columns[0]) {
        if (item.base == annos.tumor_gt) {
            normal_tumor_allele_count++;
        }
    }
    annos.cnt_tumor[0] = normal_tumor_allele_count;
    std::fill(annos.cnt_type_strand.begin(), annos.cnt_type_strand.end(), 0);
    for (int i = 0; i < n_tumor_sample; ++i) {
        bool has_mutation;
        if (is_ignore0) {
            has_mutation = (!is_empty[i] || !is_all_ref[i])
                           && (likelihoods[i][tumor_gt_idx][0] < llh_integral[i][tumor_gt_idx]);
        } else {
            has_mutation = !is_empty[i] && (likelihoods[i][tumor_gt_idx][0] < llh_integral[i][tumor_gt_idx]);
        }
        if (has_mutation) {
            annos.genotype[i + 1] = 1;
        } else {
            annos.genotype[i + 1] = 0;
        }
        annos.cnt_tumor[i + 1] = n_tumor[i * tumor_bases.size() + tumor_gt_idx];
        annos.zq[i + 1] = -10 * llh_integral[i][tumor_gt_idx];
    }
    for (int i = 0; i < n_tumor_sample + 1; ++i) {
        for (auto &read : columns[i]) {
            if (annos.normal_gt.contain(read.base)) {
                if (!read.is_reverse_strand) {
                    annos.cnt_type_strand[i*4] += 1;        // ref forward
                } else {
                    annos.cnt_type_strand[i*4 + 1] += 1;    // ref reverse
                }
            } else if (annos.tumor_gt == read.base) {
                if (!read.is_reverse_strand) {
                    annos.cnt_type_strand[i*4 + 2] += 1;    // alt forward
                } else {
                    annos.cnt_type_strand[i*4 + 3] += 1;    // alt reverse
                }
            }
        }
    }

    annos.quality = log_sum_exp_final(max_nume_elem, nume) - log_sum_exp_final(max_deno_elem, deno);
    return annos.quality;
}

double moss::qphred2prob(int qphred) {
    return pow(10.f, -static_cast<double>(qphred) / 10.f);
}

double moss::log_trinomial(unsigned long s, unsigned long k, unsigned long t) {
    double log_tri = std::lgamma(s + k + t + 1) - std::lgamma(s + 1) - std::lgamma(k + 1) - std::lgamma(t + 1);
    return log_tri;
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

template<typename T>
T moss::log_sum_exp(T a, T b) {
    if (a >= b) {
        return a + log(exp(b - a) + 1.0);
    } else {
        return b + log(exp(a - b) + 1.0);
    }
}

template<typename T>
void moss::log_sum_exp_init(T &max_elem, T &accum) {
    max_elem = -std::numeric_limits<T>::infinity();
    accum = T{};
}

template<typename T>
void moss::log_sum_exp_iter(T &max_elem, T &accum, T item) {
    if (item >= max_elem) {
        accum *= exp(max_elem - item);
        accum += 1.f;
        max_elem = item;
    } else {
        accum += exp(item - max_elem);
    }
}

template<typename T>
T moss::log_sum_exp_final(T &max_elem, T &accum) {
    return max_elem + log(accum);
}
