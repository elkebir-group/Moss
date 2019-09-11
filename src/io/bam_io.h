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

    struct _mplp_conf_t {
        int min_mq, flag, min_baseQ, capQ_thres, max_depth, max_indel_depth, fmt_flag, all;
        int rflag_require, rflag_filter;
        int openQ, extQ, tandemQ, min_support; // for indels
        double min_frac; // for indels
        char *reg, *pl_list, *fai_fname, *output_fname;
        void *bed, *rghash;
        int argc;
        char **argv;
    };
    using mplp_conf_t = _mplp_conf_t;

    struct _data_t {
        samFile *sam_fp;
        hts_itr_multi_t *iter;
        hts_idx_t *index;
        bam_hdr_t *header;
    };
    using data_t = _data_t;

    class BamStreamer {
    private:
        const unsigned long num_samples;
        std::vector<data_t **> meta;
        const int min_base_qual;
        const int min_map_qual;
        const std::string reference;
        faidx_t *ref_fp;
        const MapContigLoci loci;                           //!< candidate positions
        std::vector<std::vector<hts_reglist_t>> regions;    //!< regions of interest for sam_itr_regions
        std::vector<std::pair<locus_t, locus_t>> windows;   //!< sliding windows for each sample
        std::vector<int> tids;                              //!< current hts tid (current chrom) for each sample
        std::vector<std::set<locus_t>::iterator> iters;     //!< iterator point to loci
        std::vector<std::deque<locus_t>> actives;           //!< active loci for each sample
        std::vector<Buffer> buffers;                        //!< buffer for Pileups in building
        const static uint16_t FAIL_FLAGS = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;

    public:
        explicit BamStreamer(std::string ref_file_name,
                             const std::vector<std::string> &bam_file_names,
                             const MapContigLoci &loci,
                             int min_baseQ = 13,
                             int min_mapQ = 30);

        virtual ~BamStreamer();

        Pileups get_column();

        int print_pileup_seq(const bam_pileup1_t *p, int n);
    };
}


#endif //MOSS_BAM_IO_H
