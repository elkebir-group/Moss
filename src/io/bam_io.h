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

    struct data_t {
        samFile *sam_fp;
        hts_itr_multi_t *iter;
        hts_idx_t *index;
        bam_hdr_t *header;
    };


    class BamStreamer {
    public:
        virtual std::vector<Read> get_column(std::string &ret_contig, int &ret_pos) = 0;
        virtual ~BamStreamer() {};
    };

    class SingleBamStreamer : public BamStreamer {
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
        explicit SingleBamStreamer(const std::string &ref_file_name,
                                   const std::string &bam_file_name,
                                   const MapContigLoci loci,
                                   int min_baseQ,
                                   int min_mapQ,
                                   bool filter_edit_distance);

        ~SingleBamStreamer() override;

        std::vector<Read> get_column(std::string &ret_contig, int &ret_pos) override;
    };


    class PairedBamStreamer : public BamStreamer {
    private:
        SingleBamStreamer original;
        SingleBamStreamer realigned;
    public:
        ///
        /// @param ref_file_name
        /// @param original_bam_file_name
        /// @param realigned_bam_file_name
        /// @param loci
        /// @param min_baseQ
        /// @param min_mapQ
        /// @param filter_edit_distance
        PairedBamStreamer(const std::string &ref_file_name,
                          const std::string &original_bam_file_name,
                          const std::string &realigned_bam_file_name,
                          MapContigLoci loci,
                          int min_baseQ,
                          int min_mapQ,
                          bool filter_edit_distance);

        std::vector<Read> get_column(std::string &ret_contig, int &ret_pos) override;
    };


    class MultiBamStreamer {
    private:
        const unsigned long num_samples;
        const int min_base_qual;
        const int min_map_qual;
        const std::string reference;
        faidx_t *ref_fp;
        const MapContigLoci loci;                           //!< candidate positions
        std::vector<std::unique_ptr<BamStreamer>> streams;                   //!< single bam file
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
    };
}


#endif //MOSS_BAM_IO_H
