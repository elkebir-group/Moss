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
        const int min_base_qual;
        const int min_map_qual;
        const MapContigLoci loci;                           //!< candidate positions
        data_t **meta;
        std::vector<hts_reglist_t> region;    //!< regions of interest for sam_itr_regions
        std::pair<locus_t, locus_t> window;   //!< sliding windows for each sample
        int tid;                              //!< current hts tid (current chrom) for each sample
        std::map<locus_t, Aggregate>::const_iterator iter;  //!< iterator point to loci
        std::deque<locus_t> actives;           //!< active loci for each sample
        Buffer buffer;
        bool is_filter_edit_distance;
        const static uint16_t FAIL_FLAGS = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    public:
        explicit BamStreamer(const std::string &ref_file_name,
                             const std::string &bam_file_name,
                             const MapContigLoci loci,
                             int min_baseQ,
                             int min_mapQ,
                             bool filter_edit_distance);

        // BamStreamer(const BamStreamer &other);

        virtual ~BamStreamer();

        std::vector<Read> get_column(std::string &ret_contig, int &ret_pos);
    };


    class PairedBamStreamer {
    private:
        BamStreamer original;
        BamStreamer realigned;
    public:
        PairedBamStreamer(const std::string &ref_file_name,
                          const std::string &original_bam_file_name,
                          const std::string &realigned_bam_file_name,
                          const MapContigLoci loci,
                          int min_baseQ,
                          int min_mapQ,
                          bool filter_edit_distance);

        std::vector<Read> get_column(std::string &ret_contig, int &ret_pos);
    };


    class MultiBamStreamer {
    private:
        const unsigned long num_samples;
        const int min_base_qual;
        const int min_map_qual;
        const std::string reference;
        faidx_t *ref_fp;
        const MapContigLoci loci;                           //!< candidate positions
        std::vector<PairedBamStreamer> streams;                   //!< single bam file
        const static uint16_t FAIL_FLAGS = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
        bool is_filter_edit_distance;

    public:
        explicit MultiBamStreamer(std::string ref_file_name,
                             const std::vector<std::string> &bam_file_names,
                             const std::vector<std::string> &realigned_file_names,
                             const MapContigLoci &loci,
                             int min_baseQ = 13,
                             int min_mapQ = 30,
                             bool filter_edit_distance = false);

        virtual ~MultiBamStreamer();

        Pileups get_column();

        int print_pileup_seq(const bam_pileup1_t *p, int n);
    };
}


#endif //MOSS_BAM_IO_H
