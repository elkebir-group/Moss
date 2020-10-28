//
// Created by Chuanyi Zhang on 2019-06-27.
//

#include <utility>
#include <fstream>
#include <sstream>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "vcf_io.h"

using namespace moss;

VcfWriter::VcfWriter(const std::string &filename, const MapContigLoci &loci, unsigned long num_tumor_samples,
                     const std::string& ref_file, const std::vector<std::string> &bam_files, const std::string &command,
                     bool filter_total_dp,
                     bool filter_vaf, float qual_thr)
    : is_filter_total_dp(filter_total_dp), is_filter_vaf(filter_vaf), qual_thr(qual_thr) {
    ofile = bcf_open(filename.c_str(), "w");
    header = bcf_hdr_init("w");
    // FILTER
    addFilters()
        (true, filter_low_normal_depth,
         "##FILTER=<ID=LOW_NORMAL_DP,Description=\"The read depth of normal sample is less than 6\">")
        (true, filter_low_qual, qual_thr,
         "##FILTER=<ID=LOW_QUAL,Description=\"QUAL is below the threshold\">")
        (true, filter_low_tumor_support,
         "##FILTER=<ID=LOW_TUMOR_SUPP,Description=\"Less than 4 reads support in all tumor samples\">")
        (is_filter_total_dp, filter_low_total_depth,
         "##FILTER=<ID=LOW_TOTAL_DP,Description=\"Total depth of 23 samples < 150\">")
        (is_filter_vaf, filter_low_vaf,
         "##FILTER=<ID=LOW_VAF,Description=\"VAF in any tumor samples < 0.1\">")
        (true, filter_empty_strand,
         "##FILTER=<ID=EMPTY_STRAND,Description=\"Total somatic allele count is 0 on one strand\">");
    // FORMAT
    bcf_hdr_append(header,
                   "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this sample\">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=TCOUNT,Number=1,Type=Integer,Description=\"Number of reads with alternate allele at this position for this sample\">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise"
                   "the Fisher's Exact Test to detect strand bias.\">");
    // INFO
    bcf_hdr_append(header,
                   "##INFO=<ID=TIN,Number=1,Type=Float,Description=\"Probability of normal also has tumor genotype\">");
    bcf_hdr_append(header,
                   "##INFO=<ID=SGT,Number=1,Type=String,Description=\"Most likely somatic genotype\">");
    bcf_hdr_append(header,
                   "##INFO=<ID=NUMPASS,Number=1,Type=Integer,Description=\"Number of samples that pass the base caller\">");
    ref_idx = fai_load(ref_file.c_str());
    for (const auto &contig : loci) {
        int len = faidx_seq_len(ref_idx, contig.first.c_str());
        bcf_hdr_append(header, ("##contig=<ID=" + contig.first + ",length=" + std::to_string(len) + ">").c_str());
    }
    // Miscellanous
    bcf_hdr_append(header, ("##reference=" + ref_file).c_str());
    bcf_hdr_append(header, "##program=Moss");
    auto now = std::time(nullptr);
    std::stringstream ss;
    ss << "##datetime=" << std::put_time(std::localtime(&now), "%FT%T%z");
    bcf_hdr_append(header, ss.str().c_str());
    bcf_hdr_append(header, ("##command=" + command).c_str());
    int i = 0;
    for (const auto &bam_file : bam_files) {
        samFile *fp = sam_open(bam_file.c_str(), "r");
        bam_hdr_t *hdr = sam_hdr_read(fp);
        std::string header_text{hdr->text};
        std::string name{};
        auto read_group_pos = header_text.find("@RG");
        if (read_group_pos != std::string::npos) {
            auto sample_pos = header_text.find("SM:", read_group_pos);
            if (sample_pos != std::string::npos) {
                auto sample_end = header_text.find_first_of("\t\n", sample_pos);
                name = header_text.substr(sample_pos + 3, sample_end - sample_pos - 3);
            } else {
                auto id_pos = header_text.find("ID:", read_group_pos);
                if (id_pos != std::string::npos) {
                    auto id_end = header_text.find_first_of("\t\n", id_pos);
                    name = header_text.substr(id_pos + 3, id_end - id_pos - 3);
                }
            }
        }
        if (name.empty()) {
            bcf_hdr_add_sample(header, ("sample" + std::to_string(i)).c_str());
        } else {
            bcf_hdr_add_sample(header, name.c_str());
        }
        bam_hdr_destroy(hdr);
        sam_close(fp);

        ++i;
    }
    switch (bcf_hdr_write(ofile, header)) {
        case 0:
            break;
        case -1:
        default:
            throw std::runtime_error("Failed to write to file " + filename);
    }
    rec = bcf_init();
}

VcfWriter::~VcfWriter() {
    bcf_close(ofile);
    if (ref_idx != nullptr) {
        fai_destroy(ref_idx);
    }
    bcf_hdr_destroy(header);
    bcf_destroy(rec);
}

bool VcfWriter::filter_low_qual(Annotation &anno, float thr) {
    return anno.quality <= thr;     // not multiplied by -10 yet
}

bool VcfWriter::filter_low_normal_depth(Annotation &anno) {
    return anno.cnt_read[0] >= LEAST_NORMAL_DEPTH;
}

bool VcfWriter::filter_low_tumor_support(Annotation &anno) {
    auto max_supp = std::max_element(anno.cnt_tumor.cbegin() + 1, anno.cnt_tumor.cend());
    return *max_supp >= LEAST_TUMOR_SUPPORT;
}

bool VcfWriter::filter_low_total_depth(Annotation &anno) {
    int total_depth = 0;
    for (int i = 1; i < anno.cnt_read.size(); i++) {
        total_depth += anno.cnt_read[i];
    }
    return total_depth >= LEAST_TOTAL_DEPTH;
}

bool VcfWriter::filter_low_vaf(Annotation &anno) {
    double max_vaf = 0;
    double vaf;
    for (int i = 1; i < anno.cnt_read.size(); i++) {
        if (anno.genotype[i] == 1) {
            if (anno.cnt_read[i] != 0) {
                vaf = double(anno.cnt_tumor[i]) / double(anno.cnt_read[i]);
            } else {
                vaf = 0;
            }
            if (vaf > max_vaf)
                max_vaf = vaf;
        }
    }
    return max_vaf >= LEAST_VAF;
}

bool VcfWriter::filter_empty_strand(Annotation &anno) {
    int forward_somatic_cnt = 0;
    int reverse_somatic_cnt = 0;
    for (int i = 1; i < anno.cnt_read.size(); i++) {
        forward_somatic_cnt += anno.cnt_type_strand[4 * i + 2];
        reverse_somatic_cnt += anno.cnt_type_strand[4 * i + 3];
    }
    return (forward_somatic_cnt != 0 && reverse_somatic_cnt != 0);
}

void
VcfWriter::write_record(std::string chrom, std::pair<const locus_t, Aggregate> consensus, uint8_t ref,
                        Annotation &annos, int num_samples) {
    auto &depth = annos.cnt_read;
    auto &tumor_count = annos.cnt_tumor;
    auto &zq = annos.zq;
    auto &Z = annos.genotype;
    bcf_update_filter(header, rec, nullptr, 0);
    rec->qual = -10 * annos.quality;
    rec->pos = consensus.first;
    rec->rid = bcf_hdr_name2id(header, chrom.c_str());
    uint8_t alt = annos.tumor_gt;
    std::string alleles{seq_nt16_str[ref]};
    alleles += ",";
    alleles += alt == 0 ? '.' : seq_nt16_str[alt];
    bcf_update_alleles_str(header, rec, alleles.c_str());
    bcf_update_format_int32(header, rec, "DP", depth.data(), num_samples);
    bcf_update_format_int32(header, rec, "TCOUNT", tumor_count.data(), num_samples);
    bcf_update_format_int32(header, rec, "Z", Z.data(), num_samples);
    bcf_update_format_float(header, rec, "ZQ", zq.data(), num_samples);
    bcf_update_info_float(header, rec, "TIN", &annos.log_t_in_normal, 1);
    bcf_update_info_int32(header, rec, "NUMPASS", &consensus.second.num_pass, 1);
    std::vector<float> SOR(num_samples);
    for (int i = 0; i < num_samples; i++) {
        float a = annos.cnt_type_strand[4 * i] + 1;
        float b = annos.cnt_type_strand[4 * i + 1] + 1;
        float c = annos.cnt_type_strand[4 * i + 2] + 1;
        float d = annos.cnt_type_strand[4 * i + 3] + 1;
        SOR[i] = std::log(((a * d) / (b * c) + (b * c) / (a * d))
                          * std::min(a, b) * std::max(c, d) / std::min(c, d) / std::max(a, b));
    }
    bcf_update_format_float(header, rec, "SOR", SOR.data(), num_samples);
    bcf_update_format_int32(header, rec, "SB", annos.cnt_type_strand.data(), num_samples * 4);
    bool all_clear = true;

    for (auto &&filter : filters) {
        if (!(filter.is_filter(annos))) {
            all_clear = false;
            bcf_add_filter(header, rec, filter.filter_id);
        }
    }
    if (all_clear) {
        bcf_add_filter(header, rec, filter_pass_id);
    }

    switch (bcf_write(ofile, header, rec)) {
        case 0:
            break;
        case -1:
        default:
            throw std::runtime_error(std::string("Failed to write record to ") + ofile->fn);
    }
    bcf_clear(rec);
}

FilterHelper VcfWriter::addFilters() {
    bcf_hdr_append(header, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    filter_pass_id = bcf_hdr_id2int(header, BCF_DT_ID, "PASS");
    return FilterHelper(*this);
}
