//
// Created by Chuanyi Zhang on 2019-06-27.
//

#ifndef MOSS_VCF_IO_H
#define MOSS_VCF_IO_H

#include <string>
#include <vector>
#include "../core/types.h"

namespace moss {
    using filter_t = uint32_t;

    class Vcf {
    private:
        std::string filename;
        std::vector<std::string> chrom;
        std::vector<locus_t> pos;
        std::vector<std::string> id;
        std::vector<uint8_t> ref;
        std::vector<std::vector<uint8_t>> alt;
        std::vector<BaseSet> gt;
        std::vector<uint32_t> qual;
        std::vector<std::string> filter;
        std::vector<bool> is_snv;
    public:
        explicit Vcf(const std::string &filename);

        const std::vector<locus_t> &get_pos() const;

        const std::vector<std::string> &get_filter() const;

        const std::vector<BaseSet> &get_gt() const;

        const std::vector<std::string> &get_chrom() const;
    };
}


#endif //MOSS_VCF_IO_H
