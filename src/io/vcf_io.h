//
// Created by Chuanyi Zhang on 2019-06-27.
//

#ifndef MOSS_VCF_IO_H
#define MOSS_VCF_IO_H

#include <string>
#include <utility>
#include <vector>
#include <map>
#include <functional>
#include <htslib/hfile.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "../core/types.h"

namespace moss {
    const unsigned char seq_nt16_table[256] = {
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        1, 2, 4, 8, 15, 15, 15, 15, 15, 15, 15, 15, 15, 0 /*=*/, 15, 15,
        15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15,
        15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15,
        15, 1, 14, 2, 13, 15, 15, 4, 11, 15, 15, 12, 15, 3, 15, 15,
        15, 15, 5, 6, 8, 15, 7, 9, 15, 10, 15, 15, 15, 15, 15, 15,

        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15,
        15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
    };

    using filter_t = uint32_t;

    struct RecData {
        BaseSet bases;
        bool is_pass;

        RecData() : bases(BaseSet(0xff)), is_pass(false) {};

        RecData(BaseSet bases, bool is_pass) : bases(bases), is_pass(is_pass) {};
    };

    struct LociRecData {
        int num_pass;

        LociRecData() : num_pass(-1) {};

        explicit LociRecData(int num_pass) : num_pass(num_pass) {};
    };

    template<typename T>
    class VcfReader {
    private:
        std::string filename;
        std::map<std::string, std::map<locus_t, T>> records;

        static std::string get_mode(htsFormat *format);

    public:
        /**
         * Read a VCF file and store records.
         * @param filename is the VCF file to be read.
         */
        explicit VcfReader(const std::string &filename);

        ~VcfReader();

        /**
         * Find a record.
         * @param contig is the name of the contig.
         * @param pos is the position in the contig.
         * @return The record if found, otherwise returns an empty record.
         */
        T find(const std::string &contig, locus_t pos);

        /// Get access to the records in VCF.
        /// @return Const reference to the records.
        const std::map<std::string, std::map<locus_t, T>> &get_records() const;

        /// Whether the VCF file is empty.
        /// @return Indicator of whether the VCF file is empty.
        bool empty();
    };

    template<>
    inline VcfReader<RecData>::VcfReader(const std::string &filename) : filename(filename) {
        if (!filename.empty()) {
            hFILE *fp = hopen(filename.c_str(), "r");
            htsFormat format;
            hts_detect_format(fp, &format);
            switch (hclose(fp)) {
                case 0:
                    break;
                case EOF:
                default:
                    throw std::runtime_error("Failed to close file.");
            }
            bcf_hdr_t *header;
            htsFile *ifile;
            ifile = bcf_open(filename.c_str(), get_mode(&format).c_str());
            if (ifile != nullptr) {
                header = bcf_hdr_read(ifile);
                bcf1_t *rec = bcf_init();
                int ngt{-1},
                    *gt = nullptr,
                    ngt_arr = 0;
                uint8_t normal_gt;
                while (bcf_read(ifile, header, rec) == 0) {
                    if (bcf_is_snp(rec)) {
                        const char *contig = bcf_hdr_id2name(header, rec->rid);
                        ngt = bcf_get_format_int32(header, rec, "GT", &gt, &ngt_arr);
                        if (ngt == 1 && gt[0] != 0) {
                            normal_gt = seq_nt16_table[*rec->d.allele[0]];
                            if (*rec->d.allele[bcf_gt_allele(gt[0])] != '*') {
                                normal_gt |= seq_nt16_table[*rec->d.allele[bcf_gt_allele(gt[0])]];
                            }
                        } else if (ngt == 2 && gt[0] != 0) {
                            normal_gt = 0;
                            if (*rec->d.allele[bcf_gt_allele(gt[0])] != '*') {
                                normal_gt |= seq_nt16_table[*rec->d.allele[bcf_gt_allele(gt[0])]];
                            }
                            if (*rec->d.allele[bcf_gt_allele(gt[1])] != '*') {
                                normal_gt |= seq_nt16_table[*rec->d.allele[bcf_gt_allele(gt[1])]];
                            }
                        } else {
                            normal_gt = 0xff_8;
                        }
                        bool is_rec_pass = bcf_has_filter(header, rec, "PASS");
                        auto search = records.find(contig);
                        if (search != records.end()) {
                            search->second.emplace(std::make_pair(rec->pos, RecData{BaseSet{normal_gt}, is_rec_pass}));
                        } else {
                            records.emplace(std::make_pair(contig,
                                                           std::map<locus_t, RecData>{{rec->pos, RecData{
                                                               BaseSet{normal_gt},
                                                               is_rec_pass}}}));
                        }
                    } else {
                        continue;
                    }
                }
                if (rec != nullptr) {
                    bcf_destroy(rec);
                }
                if (ngt >= 0) {
                    free(gt);
                }
                bcf_close(ifile);
                bcf_hdr_destroy(header);
            } else {
                exit(1);
            }
        }
    }

    template<>
    inline VcfReader<LociRecData>::VcfReader(const std::string &filename) : filename(filename) {
        if (!filename.empty()) {
            hFILE *fp = hopen(filename.c_str(), "r");
            htsFormat format;
            hts_detect_format(fp, &format);
            switch (hclose(fp)) {
                case 0:
                    break;
                case EOF:
                default:
                    throw std::runtime_error("Failed to close file " + filename);
            }
            bcf_hdr_t *header;
            htsFile *ifile;
            ifile = bcf_open(filename.c_str(), get_mode(&format).c_str());
            if (ifile != nullptr) {
                header = bcf_hdr_read(ifile);
                bcf1_t *rec = bcf_init();
                while (bcf_read(ifile, header, rec) == 0) {
                    if (bcf_is_snp(rec)) {
                        const char *contig = bcf_hdr_id2name(header, rec->rid);
                        int *num_pass = nullptr;
                        int n_num_pass = 0;
                        int ret = bcf_get_info_int32(header, rec, "NUMPASS", &num_pass, &n_num_pass);
                        if (ret < 0) {
                            num_pass = new int(0);
                        }
                        auto search = records.find(contig);
                        if (search != records.end()) {
                            search->second.emplace(std::make_pair(rec->pos, LociRecData{*num_pass}));
                        } else {
                            records.emplace(std::make_pair(contig,
                                                           std::map<locus_t, LociRecData>{{rec->pos, LociRecData{
                                                               *num_pass}}}));
                        }
                        free(num_pass);
                    } else {
                        continue;
                    }
                }
                if (rec != nullptr) {
                    bcf_destroy(rec);
                }
                bcf_close(ifile);
                bcf_hdr_destroy(header);
            } else {
                exit(1);
            }
        }
    }

    template<typename T>
    VcfReader<T>::~VcfReader() = default;

    template<typename T>
    T VcfReader<T>::find(const std::string &contig, locus_t pos) {
        auto search_contig = records.find(contig);
        if (search_contig != records.end()) {
            auto search_pos = search_contig->second.find(pos);
            if (search_pos != search_contig->second.end()) {
                return search_pos->second;
            } else {
                return T();
            }
        } else {
            return T();
        }
    }

    template<typename T>
    std::string VcfReader<T>::get_mode(htsFormat *format) {
        switch (format->format) {
            case bcf:
                return "rb";

            case vcf:
                if (format->compression == gzip) {
                    return "rz";
                } else {
                    return "r";
                }

            default:
                return "r";
        }
    }

    template<typename T>
    const std::map<std::string, std::map<locus_t, T>> &VcfReader<T>::get_records() const {
        return records;
    }

    template<typename T>
    bool VcfReader<T>::empty() {
        return filename.empty();
    }

    class FilterHelper;

    class VcfWriter;

    typedef std::function<bool(Annotation &)> FilterFunc;
    typedef std::function<bool(Annotation &, float)> FilterFuncWithFloat;

    class VcfWriter {
        friend class FilterHelper;

    private:
        bcf_hdr_t *header;
        bcf1_t *rec;
        htsFile *ofile;
        std::string reference;
        faidx_t *ref_idx;
        float qual_thr;

        static bool filter_low_qual(Annotation &anno, float thr);

        static bool filter_low_normal_depth(Annotation &anno);

        static bool filter_low_tumor_support(Annotation &anno);

        static bool filter_low_total_depth(Annotation &anno);

        static bool filter_low_vaf(Annotation &anno);

        static bool filter_empty_strand(Annotation &anno);

        static const int LEAST_NORMAL_DEPTH = 6;
        static const int LEAST_TUMOR_SUPPORT = 4;
        static const int LEAST_TOTAL_DEPTH = 150;
        static constexpr float LEAST_VAF = 0.1;
        int filter_pass_id{};
        bool is_filter_total_dp;
        bool is_filter_vaf;
        struct Filter {
            FilterFunc is_filter;
            int filter_id;
        };
        std::vector<Filter> filters;

        FilterHelper addFilters();

    public:
        /**
         * Construct a writer of VCF file.
         * @param filename is the name of VCF output.
         * @param loci contains the contigs.
         * @param num_tumor_samples is the number of tumor samples.
         * @param ref_file is the name of the reference FASTA file.
         * @param bam_files are a vector of BAM files.
         * @param command is the concatenated command of Moss.
         * @param filter_total_dp indicates whether filter on total depth.
         * @param filter_vaf indicates whether filter on minimum VAF.
         * @param qual_thr is the threshold of quality score, default is 0 (unused);
         */
        VcfWriter(const std::string &filename, const MapContigLoci &loci, unsigned long num_tumor_samples,
                  const std::string &ref_file, const std::vector<std::string> &bam_files, const std::string &command,
                  bool filter_total_dp = false,
                  bool filter_vaf = false, float qual_thr = 0);

        ~VcfWriter();

        /**
         * Write a record to the output VCF file.
         *
         * @param chrom is the record's chromosome.
         * @param consensus is the aggregated information.
         * @param ref is the reference base.
         * @param annos is the annotations.
         * @param num_samples is the total number of samples, including normal.
         */
        void
        write_record(std::string chrom, std::pair<const locus_t, Aggregate> consensus, uint8_t ref, Annotation &annos,
                     int num_samples);
    };

    class FilterHelper {
    public:
        FilterHelper(VcfWriter &writer) : writer(writer) {}

        FilterHelper &operator()(bool toggle, FilterFunc func, const std::string &line) {
            if (toggle) {
                bcf_hdr_append(writer.header, line.c_str());
                auto id_start_pos = line.find("ID") + 3;
                auto id_end_pos = line.find_first_of(',', id_start_pos);
                std::string id = line.substr(id_start_pos, id_end_pos - id_start_pos);
                int filter_id = bcf_hdr_id2int(writer.header, BCF_DT_ID, id.c_str());
                writer.filters.push_back({std::move(func), filter_id});
            }
            return *this;
        }

        FilterHelper &operator()(bool toggle, const FilterFuncWithFloat &func, float thr, const std::string &line) {
            if (toggle) {
                bcf_hdr_append(writer.header, line.c_str());
                auto id_start_pos = line.find("ID") + 3;
                auto id_end_pos = line.find_first_of(',', id_start_pos);
                std::string id = line.substr(id_start_pos, id_end_pos - id_start_pos);
                int filter_id = bcf_hdr_id2int(writer.header, BCF_DT_ID, id.c_str());
                writer.filters.push_back({std::bind(func, std::placeholders::_1, thr), filter_id});
            }
            return *this;
        }

    private:
        VcfWriter &writer;
    };
}


#endif //MOSS_VCF_IO_H
