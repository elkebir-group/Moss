//
// Created by Chuanyi Zhang on 2019-06-27.
//

#include <utility>
#include <fstream>
#include <sstream>
#include <htslib/hfile.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <algorithm>
#include "vcf_io.h"

using namespace moss;

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

VcfReader::VcfReader(const std::string &filename) : filename(filename) {
    if (!filename.empty()) {
        hFILE *fp = hopen(filename.c_str(), "r");
        htsFormat format;
        hts_detect_format(fp, &format);
        hclose(fp);
        bcf_hdr_t *header;
        bcf1_t *rec;
        htsFile *ifile;
        ifile = bcf_open(filename.c_str(), get_mode(&format).c_str());
        if (ifile != nullptr) {
            header = bcf_hdr_read(ifile);
            bcf1_t *rec = bcf_init();
            int ngt,
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
                                                       std::map<locus_t, RecData>{{rec->pos, RecData{BaseSet{normal_gt},
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

VcfReader::~VcfReader() = default;

RecData VcfReader::find(const std::string &contig, locus_t pos) {
    auto search_contig = records.find(contig);
    if (search_contig != records.end()) {
        auto search_pos = search_contig->second.find(pos);
        if (search_pos != search_contig->second.end()) {
            return search_pos->second;
        } else {
            return RecData{BaseSet(0xff), false};
        }
    } else {
        return RecData{BaseSet(0xff), false};
    }
}

std::string VcfReader::get_mode(htsFormat *format) {
    switch (format->format) {
        case bcf:
            return "rb";
            break;

        case vcf:
            if (format->compression == gzip) {
                return "rz";
            } else if (format->compression == no_compression) {
                return "r";
            } else {
                return "r";
            }
            break;

        default:
            return "r";
            break;
    }
}

const std::map<std::string, std::map<locus_t, RecData>> &VcfReader::get_records() const {
    return records;
}

bool VcfReader::empty() {
    return filename.empty();
}

VcfWriter::VcfWriter(const std::string &filename, MapContigLoci loci, unsigned long num_tumor_samples,
                     std::string ref_file, std::vector<std::string> bam_files, bool filter_total_dp, bool filter_vaf, float qual_thr)
    : is_filter_total_dp(filter_total_dp), qual_thr(qual_thr) {
    ofile = bcf_open(filename.c_str(), "w");
    header = bcf_hdr_init("w");
    // FILTER
    addFilters()
    (true, filter_low_normal_depth,
        "##FILTER=<ID=LOW_NORMAL_DP,Description=\"The read depth of normal sample is less than 6\">")
    (true, filter_low_qual, qual_thr,
        "##FILTER=<ID=LOW_QUAL,Description=\"QUAL is below the threshold\">")
    (true, filter_low_tumor_support,
        "##FILTER=<ID=LOW_TUMOR_SUPP,Description=\"Less than 4 reads support any tumor sample\">")
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
                   "##FORMAT=<ID=Z,Number=1,Type=Integer,Description=\"Boolean indicator of containing somatic SNV in this sample \">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=ZQ,Number=1,Type=Float,Description=\"Sample-specific likelihood of boolean indicator to be true (phred-scaled)\">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=SB,Number=4,Type=Integer,Description=\"Per-sample component statistics which comprise"
                   "the Fisher's Exact Test to detect strand bias.\">");
    // INFO
    bcf_hdr_append(header,
                   "##INFO=<ID=TIN,Number=1,Type=Float,Description=\"Probability of normal also has tumor genotype\">");
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
    bcf_hdr_append(header, "##command=");
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
                auto sample_end = header_text.find('\t', sample_pos);
                name = header_text.substr(sample_pos + 3, sample_end - sample_pos - 3);
            } else {
                auto id_pos = header_text.find("ID:", read_group_pos);
                if (id_pos != std::string::npos) {
                    auto id_end = header_text.find('\t', id_pos);
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
    bcf_hdr_write(ofile, header);
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
        forward_somatic_cnt += anno.cnt_type_strand[4*i + 2];
        reverse_somatic_cnt += anno.cnt_type_strand[4*i + 3];
    }
    return (forward_somatic_cnt != 0 && reverse_somatic_cnt != 0);
}

void
VcfWriter::write_record(std::string chrom, int pos, uint8_t ref, Annotation &annos, int num_samples) {
    auto &depth = annos.cnt_read;
    auto &tumor_count = annos.cnt_tumor;
    auto &zq = annos.zq;
    auto &Z = annos.genotype;
    bcf_update_filter(header, rec, nullptr, 0);
    rec->qual = -10 * annos.quality;
    rec->pos = pos;
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
    std::vector<float> SOR(num_samples);
    for (int i = 0; i < num_samples; i++) {
        float a = annos.cnt_type_strand[4*i] + 1;
        float b = annos.cnt_type_strand[4*i + 1] + 1;
        float c = annos.cnt_type_strand[4*i + 2] + 1;
        float d = annos.cnt_type_strand[4*i + 3] + 1;
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
    
    bcf_write(ofile, header, rec);
    bcf_clear(rec);
}

FilterHelper VcfWriter::addFilters() {
    bcf_hdr_append(header, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    filter_pass_id = bcf_hdr_id2int(header, BCF_DT_ID, "PASS");
    return FilterHelper(*this);
}
