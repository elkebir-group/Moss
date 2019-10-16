//
// Created by Chuanyi Zhang on 2019-06-27.
//

#ifndef MOSS_VCF_IO_H
#define MOSS_VCF_IO_H

#include <string>
#include <vector>
#include <map>
#include "htslib/vcf.h"
#include <htslib/faidx.h>
#include "../core/types.h"

namespace moss {
    using filter_t = uint32_t;

    struct _RecData {
        BaseSet bases;
        bool is_pass;
    };

    using RecData = _RecData;
    using filter_func = bool (*)(Annotation);

    class VcfReader {
    private:
        std::string filename;
        std::map<std::string, std::map<locus_t, RecData>> records;

        static std::string get_mode(htsFormat *format);

    public:
        explicit VcfReader(const std::string &filename);

        ~VcfReader();

        RecData find(const std::string &contig, locus_t pos);

        const std::map<std::string, std::map<locus_t, RecData>> &get_records() const;

        bool empty();
    };

    class VcfWriter {
    private:
        bcf_hdr_t *header;
        bcf1_t *rec;
        htsFile *ofile;
        std::string reference;
        faidx_t *ref_idx;

        int filter_pass_id;
        int filter_low_id;
        int filter_low_normal_id;
        int filter_tumor_supp_id;
        int filter_total_dp_id;
        int filter_vaf_id;
        bool is_filter_total_dp;

    public:
        VcfWriter(const std::string &filename, MapContigLoci loci, unsigned long num_tumor_samples,
                  std::string ref_file, std::vector<std::string> bam_files, bool filter_total_dp = false);

        ~VcfWriter();

        void write_record(std::string chrom, int pos, uint8_t ref, uint8_t alt, float qual, Annotation annos, float thr, int num_samples);
    };
}


#endif //MOSS_VCF_IO_H
