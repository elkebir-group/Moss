//
// Created by Chuanyi Zhang on 2019-05-16.
//

#ifndef MOSS_BAM_IO_H
#define MOSS_BAM_IO_H

#include <vector>
#include "htslib/sam.h"
#include "../core/types.h"

namespace moss {
    using data_t = struct {
        samFile *sam_fp;
        hts_itr_t *iter;
        hts_idx_t *index;
        bam_hdr_t *header;
//    mplp_ref_t *ref;
//    const mplp_conf_t *conf;
    };


    class BamStreamer {
    private:
        unsigned long num_samples;
        std::vector<data_t **> meta;
        std::vector<bam_mplp_t> bam_handler;

    public:
        BamStreamer() = default;

        explicit BamStreamer(std::vector<std::string> bam_file_names);

        virtual ~BamStreamer();

        ReadColumns get_column(std::string contig, int locus);

        int print_pileup_seq(const bam_pileup1_t *p, int n);
    };
}


#endif //MOSS_BAM_IO_H
