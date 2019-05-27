//
// Created by Chuanyi Zhang on 2019-05-16.
//

#ifndef MOSS_BAM_IO_H
#define MOSS_BAM_IO_H

#include <vector>
#include <string>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "../core/types.h"

namespace moss {
    using plp_fp = int (*)(void *data, bam1_t *b);

    using mplp_conf_t = struct {
        int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all;
        int rflag_require, rflag_filter;
        int openQ, extQ, tandemQ, min_support; // for indels
        double min_frac; // for indels
        char *reg, *pl_list, *fai_fname, *output_fname;
//        faidx_t *fai;
        void *bed, *rghash;
        int argc;
        char **argv;
//        sam_global_args ga;
    };

    using data_t = struct {
        samFile *sam_fp;
        hts_itr_t *iter;
        hts_idx_t *index;
        bam_hdr_t *header;
//    mplp_ref_t *ref;
        const mplp_conf_t *conf;
    };

    enum class Stepper {
        nofilter,
        samtools
    };


    class BamStreamer {
    private:
        unsigned long num_samples;
        std::vector<data_t **> meta;
        std::vector<bam_mplp_t> bam_handler;
        plp_fp pileup_func;
        int min_base_qual;
        std::string reference;
        faidx_t *ref_fp;


        static int mplp_func_plain(void *data, bam1_t *b);

        static int mplp_func_samtools(void *data, bam1_t *b);

    public:
        BamStreamer() = default;

        explicit BamStreamer(std::string ref_file_name,
                             std::vector<std::string> bam_file_names,
                             Stepper stepper = Stepper::samtools,
                             int min_baseQ = 13);

        virtual ~BamStreamer();

        Pileups get_column(std::string contig, int locus);

        int print_pileup_seq(const bam_pileup1_t *p, int n);
    };
}


#endif //MOSS_BAM_IO_H
