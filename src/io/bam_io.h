//
// Created by Chuanyi Zhang on 2019-05-16.
//

#ifndef MOSS_BAM_IO_H
#define MOSS_BAM_IO_H

#include <vector>
#include <string>
#include <deque>
#include <iterator>
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "../core/types.h"

namespace moss {
    using plp_fp = int (*)(void *data, bam1_t *b);

    using Buffer = std::deque<std::vector<Read>>;

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
        hts_itr_multi_t *iter;
        hts_idx_t *index;
        bam_hdr_t *header;
//    mplp_ref_t *ref;
        // const mplp_conf_t *conf;
    };


    class BamStreamer {
    private:
        const unsigned long num_samples;
        std::vector<data_t **> meta;
        const int min_base_qual;
        const std::string reference;
        faidx_t *ref_fp;
        const MapContigLoci loci;                           //!< candidate positions
        std::vector<std::vector<hts_reglist_t>> regions;    //!< regions of interest for sam_itr_regions
        std::vector<std::pair<locus_t, locus_t>> windows;   //!< sliding windows for each sample
        std::vector<int> tids;                              //!< current hts tid (current chrom) for each sample
        std::vector<std::set<locus_t>::iterator> iters;     //!< iterator point to loci
        std::vector<std::deque<locus_t>> actives;           //!< active loci for each sample
        std::vector<Buffer> buffers;                        //!< buffer for Pileups in building

    public:
        explicit BamStreamer(std::string ref_file_name,
                             const std::vector<std::string> &bam_file_names,
                             const MapContigLoci &loci,
                             int min_baseQ = 13);

        virtual ~BamStreamer();

        Pileups get_column();

        int print_pileup_seq(const bam_pileup1_t *p, int n);
    };
}


#endif //MOSS_BAM_IO_H
