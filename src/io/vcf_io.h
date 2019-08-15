//
// Created by Chuanyi Zhang on 2019-06-27.
//

#ifndef MOSS_VCF_IO_H
#define MOSS_VCF_IO_H

#include <string>
#include <vector>
#include <map>
#include "htslib/vcf.h"
#include "../core/types.h"

namespace moss {
    using filter_t = uint32_t;

    struct _RecData {
        BaseSet bases;
        bool is_pass;
    };

    using RecData = _RecData;

    class VcfReader {
    private:
        std::string filename;
        std::map<std::string, std::map<locus_t, RecData>> records;

        std::string get_mode(htsFormat *format);

    public:
        explicit VcfReader(const std::string &filename);

        ~VcfReader();

        RecData find(const std::string &contig, locus_t pos);

        const std::map<std::string, std::map<locus_t, RecData>> &get_records() const;

        bool empty(void);
    };

    class VcfWriter {
    private:
        bcf_hdr_t *header;
        bcf1_t *rec;
        htsFile *ofile;

        int filter_pass_id;
        int filter_low_id;
    public:
        VcfWriter(const std::string &filename, MapContigLoci loci, unsigned long num_tumor_samples);

        ~VcfWriter();

        void
        write_record(std::string chrom, int pos, uint8_t ref, uint8_t alt, float qual, int *depth, int *tumor_count,
                     float thr, int num_tumor_samples);
    };
}


#endif //MOSS_VCF_IO_H
