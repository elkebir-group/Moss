//
// Created by Chuanyi Zhang on 2019-06-27.
//

#include <utility>
#include <fstream>
#include <sstream>
#include <htslib/hts.h>
#include "vcf_io.h"

using namespace moss;

const unsigned char seq_nt16_table[256] = {
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0 /*=*/,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
    15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
    15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,

    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
    15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

Vcf::Vcf(const std::string &filename) : filename(filename) {
    std::ifstream ifile{filename};
    if (ifile) {
        bool has_GT = false;
        uint32_t cnt = 0;
        while (!ifile.eof()) {
            std::string line;
            std::getline(ifile, line);
            if (line.size() > 0) {
                if (line.at(0) == '#') {
                    if (line.compare(2, 6, "FORMAT") == 0){
                        if (line.compare(13, 2, "GT") == 0) {
                            has_GT = true;
                        }
                    }
                } else {
                    std::stringstream ss{line};
                    std::string chrom_s;
                    std::string temp;
                    std::string id_s;
                    std::string ref_s;
                    std::string alt_s;
                    std::string filter_s;
                    bool snv_b = true;
                    std::getline(ss, chrom_s, '\t');
                    chrom.emplace_back(chrom_s);
                    std::getline(ss, temp, '\t');
                    pos[std::stoi(temp) - 1] = cnt;
                    cnt++;
                    std::getline(ss, id_s, '\t');
                    id.emplace_back(id_s);
                    std::getline(ss, ref_s, '\t');
                    if (ref_s.size() > 1) {
                        snv_b = false;
                    }
                    ref.emplace_back(seq_nt16_table[ref_s[0]]);
                    std::getline(ss, alt_s, '\t');
                    if (alt_s.size() > 0) {
                        std::stringstream alt_ss {alt_s};
                        std::string allele;
                        alt.emplace_back(std::vector<uint8_t>());
                        while (std::getline(alt_ss, allele, ',')) {
                            if (allele.size() == 1) {
                                alt.back().emplace_back(seq_nt16_table[allele[0]]);
                            } else {
                                snv_b = false;
                            }
                        }
                    }
                    is_snv.push_back(snv_b);
                    std::getline(ss, temp, '\t');
                    if (temp == ".") {
                        qual.emplace_back(-1);
                    } else {
                        qual.emplace_back(std::stoi(temp));
                    }
                    std::getline(ss, filter_s, '\t');
                    filter.emplace_back(filter_s);
                    std::getline(ss, temp, '\t');
                    std::getline(ss, temp, '\t');
                    if (has_GT) {
                        int idx_gt = 0;
                        {
                            std::stringstream sample{temp};
                            std::string field;
                            while (std::getline(sample, field, ':')) {
                                if (field.compare(0, 2, std::string{"GT"}) == 0) {
                                    break;
                                }
                                idx_gt++;
                            }
                        }
                        std::getline(ss, temp, '\t');
                        {
                            std::stringstream sample{temp};
                            std::string field;
                            int idx = 0;
                            while (std::getline(sample, field, ':')) {
                                if (idx == idx_gt) {
                                    if (field.size() == 3) {
                                        if ((field[0] == '0' && field[2] == '1') ||
                                            (field[0] == '1' && field[2] == '0')) {
                                            gt.emplace_back(BaseSet(ref.back() | alt.back()[0]));
                                        } else if ((field[0] == '0' && field[2] == '2') ||
                                                   (field[0] == '2' && field[2] == '0')) {
                                            gt.emplace_back(BaseSet(ref.back() | alt.back()[2]));
                                        } else if ((field[0] == '1' && field[2] == '2') ||
                                                   (field[0] == '2' && field[2] == '1')) {
                                            gt.emplace_back(BaseSet(alt.back()[0] | alt.back()[1]));
                                        } else if (field[0] == '1' && field[2] == '1') {
                                            gt.emplace_back(BaseSet(alt.back()[0]));
                                        }
                                    }
                                    break;
                                }
                                idx++;
                            }
                        }
                    }
                }

            }
        }
    }
}

const std::map<locus_t, uint32_t> &Vcf::get_pos() const {
    return pos;
}

const std::vector<std::string> &Vcf::get_filter() const {
    return filter;
}

const std::vector<BaseSet> &Vcf::get_gt() const {
    return gt;
}

const std::vector<std::string> &Vcf::get_chrom() const {
    return chrom;
}
