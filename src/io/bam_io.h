//
// Created by Chuanyi Zhang on 2019-05-16.
//

#ifndef MOSS_BAM_IO_H
#define MOSS_BAM_IO_H

#include <vector>
#include <string>
#include <deque>
#include <memory>
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
        /**
         * Build a streamer for a BAM files.
         *
         * @param ref_file_name is the name of the reference FASTA file.
         * @param bam_file_name is the name of the BAM file.
         * @param loci is the map of contig:loci to be iterated.
         * @param min_baseQ is the minimum base quality score.
         * @param min_mapQ is the minimum mapping quality score.
         * @param filter_edit_distance indicates whether filter reads whose edit distance > 3.
         */
        explicit SingleBamStreamer(const std::string &ref_file_name,
                                   const std::string &bam_file_name,
                                   MapContigLoci loci,
                                   int min_baseQ,
                                   int min_mapQ,
                                   bool filter_edit_distance);

        ~SingleBamStreamer() override;

        /**
         * Get the column of reads for the next locus in SingleBamStreamer::loci.
         * @param ret_contig is the reference of the contig string for the returned column.
         * @param ret_pos is the reference of the position for the returned column.
         * @return A vector of reads for the next locus.
         */
        std::vector<Read> get_column(std::string &ret_contig, int &ret_pos) override;
    };


    class PairedBamStreamer : public BamStreamer {
    private:
        SingleBamStreamer original;
        SingleBamStreamer realigned;
    public:
        /**
         * Build a streamer for a pair of original-realigned BAM files.
         *
         * @param ref_file_name is the name of the reference FASTA file.
         * @param original_bam_file_name is the name of the original BAM file.
         * @param realigned_bam_file_name is the name of the realigned BAM file.
         * @param loci is the map of contig:loci to be iterated.
         * @param min_baseQ is the minimum base quality score.
         * @param min_mapQ is the minimum mapping quality score.
         * @param filter_edit_distance indicates whether filter reads whose edit distance > 3.
         */
        PairedBamStreamer(const std::string &ref_file_name,
                          const std::string &original_bam_file_name,
                          const std::string &realigned_bam_file_name,
                          const MapContigLoci& loci,
                          int min_baseQ,
                          int min_mapQ,
                          bool filter_edit_distance);

        /**
         * @brief Get the column of reads for the next locus in PairedBamStreamer::loci.
         *
         * The realigned reads in the original BAM is replaced, identified by id.
         *
         * @param ret_contig is the reference of the contig string for the returned column.
         * @param ret_pos is the reference of the position for the returned column.
         * @return A vector of reads for the next locus.
         */
        std::vector<Read> get_column(std::string &ret_contig, int &ret_pos) override;
    };


    class MultiBamStreamer {
    private:
        const unsigned long num_samples;
        const int min_base_qual;
        const int min_map_qual;
        const std::string reference;
        faidx_t *ref_fp;
        const MapContigLoci loci;
        std::vector<std::unique_ptr<BamStreamer>> streams;       //!< Abstract base class of streamers
        bool is_filter_edit_distance;

    public:
        /**
         * @brief Build a streamer for multiple BAM files, can be single or pairs of original-realigned BAM files.
         *
         * `bam_file_names` and `realigned_bam_file_name` must have same length.
         *  If a BAM file is not paired, the corresponding element in `realigned_bam_file_name`
         *  will be an empty string.
         *
         * @param ref_file_name is the name of the reference FASTA file.
         * @param bam_file_names is the name of the original BAM file.
         * @param realigned_bam_file_name is the name of the realigned BAM file, can be ommitted.
         * @param loci is the map of contig:loci to be iterated.
         * @param min_baseQ is the minimum base quality score.
         * @param min_mapQ is the minimum mapping quality score.
         * @param filter_edit_distance indicates whether filter reads whose edit distance > 3.
         */
        explicit MultiBamStreamer(std::string ref_file_name,
                                  const std::vector<std::string> &bam_file_names,
                                  const std::vector<std::string> &realigned_file_names,
                                  const MapContigLoci &loci,
                                  int min_baseQ = 13,
                                  int min_mapQ = 30,
                                  bool filter_edit_distance = false);

        virtual ~MultiBamStreamer();

        /**
         * Get the column of reads for the next locus in MultiBamStreamer::loci.
         *
         * @return A vector of reads for the next locus.
         */
        Pileups get_column();
    };
}


#endif //MOSS_BAM_IO_H
