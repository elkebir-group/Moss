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
                     std::string ref_file, std::vector<std::string> bam_files, bool filter_total_dp)
    : is_filter_total_dp(filter_total_dp) {
    ofile = bcf_open(filename.c_str(), "w");
    header = bcf_hdr_init("w");
    // FILTER
    bcf_hdr_append(header, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    bcf_hdr_append(header, "##FILTER=<ID=LOW_QUAL,Description=\"QUAL is below the threshold\">");
    bcf_hdr_append(header,
                   "##FILTER=<ID=LOW_NORMAL_DP,Description=\"The read depth of normal sample is less than 6\">");
    bcf_hdr_append(header, "##FILTER=<ID=LOW_TUMOR_SUPP,Description=\"Less than 4 reads support any tumor sample\">");
    if (is_filter_total_dp) {
        bcf_hdr_append(header, "##FILTER=<ID=LOW_TOTAL_DP,Description=\"Total depth of 23 samples < 150\">");
        filter_total_dp_id = bcf_hdr_id2int(header, BCF_DT_ID, "LOW_TOTAL_DP");
    }
    bcf_hdr_append(header, "##FILTER=<ID=LOW_VAF,Description=\"VAF in any tumor samples < 0.1\">");
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
    filter_pass_id = bcf_hdr_id2int(header, BCF_DT_ID, "PASS");
    filter_low_id = bcf_hdr_id2int(header, BCF_DT_ID, "LOW_QUAL");
    filter_low_normal_id = bcf_hdr_id2int(header, BCF_DT_ID, "LOW_NORMAL_DP");
    filter_tumor_supp_id = bcf_hdr_id2int(header, BCF_DT_ID, "LOW_TUMOR_SUPP");
    filter_vaf_id = bcf_hdr_id2int(header, BCF_DT_ID, "LOW_VAF");
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

void
VcfWriter::write_record(std::string chrom, int pos, uint8_t ref, uint8_t alt, float qual, Annotation annos,
                        float thr, int num_samples) {
    auto &depth = annos.cnt_read;
    auto &tumor_count = annos.cnt_tumor;
    auto &zq = annos.zq;
    auto &Z = annos.genotype;
    bcf_update_filter(header, rec, nullptr, 0);
    rec->qual = qual;
    rec->pos = pos;
    rec->rid = bcf_hdr_name2id(header, chrom.c_str());
    std::string alleles{seq_nt16_str[ref]};
    alleles += ",";
    alleles += alt == 0 ? '.' : seq_nt16_str[alt];
    bcf_update_alleles_str(header, rec, alleles.c_str());
    bcf_update_format_int32(header, rec, "DP", depth.data(), num_samples);
    bcf_update_format_int32(header, rec, "TCOUNT", tumor_count.data(), num_samples);
    bcf_update_format_int32(header, rec, "Z", Z.data(), num_samples);
    bcf_update_format_float(header, rec, "ZQ", zq.data(), num_samples);
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
    if (depth[0] < 6) {
        bcf_add_filter(header, rec, filter_low_normal_id);
        all_clear = false;
    }
    if (qual <= thr) {
        bcf_add_filter(header, rec, filter_low_id);
        all_clear = false;
    }
    int total_cnt = 0;
    double max_vaf = 0;
    int max_supp = 0;
    for (int i = 1; i < num_samples; i++) {
        total_cnt += depth[i];
        if (Z[i] == 1) {
            double vaf;
            if (depth[i] != 0) {
                vaf = double(tumor_count[i]) / double(depth[i]);
            } else {
                vaf = 0;
            }
            if (vaf > max_vaf) {
                max_vaf = vaf;
            }
            if (tumor_count[i] > max_supp) {
                max_supp = tumor_count[i];
            }
        }
    }
    if (max_supp < 4) {
        bcf_add_filter(header, rec, filter_tumor_supp_id);
        all_clear = false;
    }
    if (is_filter_total_dp) {    
        if (total_cnt < 150) {
            bcf_add_filter(header, rec, filter_total_dp_id);
            all_clear = false;
        }
    }
    if (max_vaf < 0.1) {
        bcf_add_filter(header, rec, filter_vaf_id);
        all_clear = false;
    }
    if (all_clear) {
        bcf_add_filter(header, rec, filter_pass_id);
    }
    bcf_write(ofile, header, rec);
    bcf_clear(rec);
}
