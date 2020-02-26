//
// Created by Chuanyi Zhang on 2019-06-27.
//

#ifndef MOSS_VCF_IO_H
#define MOSS_VCF_IO_H

#include <string>
#include <vector>
#include <map>
#include <functional>
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

    class FilterHelper;
    class VcfWriter;
    typedef std::function<bool(Annotation&)> FilterFunc;
    typedef std::function<bool(Annotation&, float)> FilterFuncWithFloat;

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
        int filter_pass_id;
        bool is_filter_total_dp;
        bool is_filter_vaf;
        // std::vector<std::pair<FilterFunc, int>> filters;
        struct Filter {FilterFunc is_filter; int filter_id;};
        std::vector<Filter> filters;

    public:
        VcfWriter(const std::string &filename, const MapContigLoci &loci, unsigned long num_tumor_samples,
                  std::string ref_file, const std::vector<std::string> &bam_files, bool filter_total_dp = false, bool filter_vaf = false, float qual_thr = 0);

        ~VcfWriter();

        void write_record(std::string chrom, int pos, uint8_t ref, Annotation &annos, int num_samples);

        FilterHelper addFilters();
    };

    class FilterHelper {
    public:
        FilterHelper(VcfWriter& writer) : writer(writer) {}

        FilterHelper& operator()(bool toggle, FilterFunc func, const std::string& line) {
            if (toggle) {
                bcf_hdr_append(writer.header, line.c_str());
                auto id_start_pos = line.find("ID") + 3;
                auto id_end_pos = line.find_first_of(',', id_start_pos);
                std::string id = line.substr(id_start_pos, id_end_pos - id_start_pos);
                int filter_id = bcf_hdr_id2int(writer.header, BCF_DT_ID, id.c_str());
                writer.filters.push_back({func, filter_id});
            }
            return *this;
        }

        FilterHelper& operator()(bool toggle, FilterFuncWithFloat func, float thr, const std::string& line) {
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
