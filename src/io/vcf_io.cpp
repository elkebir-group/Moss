//
// Created by Chuanyi Zhang on 2019-06-27.
//

#include <utility>
#include <fstream>
#include <sstream>
#include <htslib/hfile.h>
#include <htslib/hts.h>
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
                    normal_gt = seq_nt16_table[*rec->d.allele[0]] |
                                seq_nt16_table[*rec->d.allele[bcf_gt_allele(gt[0])]];
                } else if (ngt == 2 && gt[0] != 0) {
                    normal_gt = seq_nt16_table[*rec->d.allele[bcf_gt_allele(gt[0])]] |
                                seq_nt16_table[*rec->d.allele[bcf_gt_allele(gt[1])]];
                } else {
                    normal_gt = 0xff_8;
                }
                bool is_rec_pass = bcf_has_filter(header, rec, "PASS");
                auto search = records.find(contig);
                if (search != records.end()) {
                    search->second.emplace(std::make_pair(rec->pos, RecData{BaseSet{normal_gt}, is_rec_pass}));
                } else {
                    records.emplace(std::make_pair(contig,
                                                   std::map<locus_t, RecData>{{rec->pos,
                                                                               RecData{BaseSet{normal_gt},
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

VcfReader::~VcfReader() {
}

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

VcfWriter::VcfWriter(const std::string &filename, MapContigLoci loci, unsigned long num_tumor_samples) {
    ofile = bcf_open(filename.c_str(), "w");
    header = bcf_hdr_init("w");
    // FILTER
    bcf_hdr_append(header, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    bcf_hdr_append(header, "##FILTER=<ID=LOW_QUAL,Description=\"QUAL is below the threshold\">");
    // FORMAT
    bcf_hdr_append(header,
                   "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at this position for this sample\">");
    bcf_hdr_append(header,
                   "##FORMAT=<ID=TCOUNT,Number=1,Type=Integer,Description=\"Number of tumor reads at this position for this sample\">");
    // TODO: get length from reference or VCF
    for (const auto &chrom : loci) {
        bcf_hdr_append(header, ("##contig=<ID=" + chrom.first + ",length=0>").c_str());
    }
    // TODO: use sample name if they have one
    for (size_t i = 0; i < num_tumor_samples; i++) {
        bcf_hdr_add_sample(header, ("sample" + std::to_string(i)).c_str());
    }
    bcf_hdr_write(ofile, header);
    filter_pass_id = bcf_hdr_id2int(header, BCF_DT_ID, "PASS");
    filter_low_id = bcf_hdr_id2int(header, BCF_DT_ID, "LOW_QUAL");
    rec = bcf_init();
}

VcfWriter::~VcfWriter() {
    bcf_hdr_destroy(header);
    bcf_destroy(rec);
    bcf_close(ofile);
}

void
VcfWriter::write_record(std::string chrom, int pos, uint8_t ref, uint8_t alt, float qual, int *depth, int *tumor_count,
                        float thr, int num_tumor_samples) {
    bcf_update_filter(header, rec, NULL, 0);
    rec->qual = qual;
    rec->pos = pos;
    rec->rid = bcf_hdr_name2id(header, chrom.c_str());
    std::string alleles{seq_nt16_str[ref]};
    alleles += ",";
    alleles += seq_nt16_str[alt];
    bcf_update_alleles_str(header, rec, alleles.c_str());
    bcf_update_format_int32(header, rec, "DP", depth, num_tumor_samples);
    bcf_update_format_int32(header, rec, "TCOUNT", tumor_count, num_tumor_samples);
    if (qual < thr) {
        bcf_update_filter(header, rec, &filter_pass_id, 1);
    } else {
        bcf_update_filter(header, rec, &filter_low_id, 1);
    }
    bcf_write(ofile, header, rec);
    bcf_clear(rec);
}
