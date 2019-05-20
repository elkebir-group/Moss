//
// Created by Chuanyi Zhang on 2019-05-16.
//

#include "bam_io.h"
#include <iostream>
#include <climits>
#include <iomanip>

using namespace moss;

BamStreamer::BamStreamer(std::vector<std::string> bam_file_names) {
    num_samples = bam_file_names.size();
//    meta = new data_t*[num_samples];
    meta.reserve(num_samples);
    bam_handler.reserve(num_samples);
    htsFormat fmt = {};
    for (const std::string &bamFile : bam_file_names) {
        samFile *fp = sam_open_format(bamFile.c_str(), "rb", &fmt);
        if (fp == nullptr) {
            std::cerr << "Error: Failed to open BAM file " << bamFile << std::endl;
            break;
        }
        hts_idx_t *idx = sam_index_load(fp, bamFile.c_str());
        if (idx == nullptr) {
            std::cerr << "Error: Failed to open BAM index " << bamFile << std::endl;
            break;
        }
        bam_hdr_t *hdr = sam_hdr_read(fp);
        if (hdr == nullptr) {
            std::cerr << "Error: Failed to fetch header from BAM file " << bamFile
                      << ".fai , the file may be broken."
                      << std::endl;
            break;
        }
        auto tmp = new data_t *[1];
        tmp[0] = new data_t[1];
        tmp[0][0] = {fp, nullptr, idx, hdr};
        meta.emplace_back(tmp);
    }
}

BamStreamer::~BamStreamer() {
    for (auto &item : meta) {
        sam_close(item[0]->sam_fp);
        bam_hdr_destroy(item[0]->header);
        if (item[0]->iter) {
            hts_itr_destroy(item[0]->iter);
        }
        delete[](item[0]);
        delete[](item);
    }
    for (auto &handler : bam_handler) {
        bam_mplp_destroy(handler);
    }
}

static int mplp_func_plain(void *data, bam1_t *b) {
    int ret;
    auto d = static_cast<data_t *>(data);
    ret = sam_itr_next(d->sam_fp, d->iter, b);
    return ret;
}

ReadColumns BamStreamer::get_column(std::string contig, int locus) {
    for (int j = 0; j < num_samples; ++j) {
        int tid = bam_name2id(meta[j][0]->header, contig.c_str());
        if ((meta[j][0]->iter = sam_itr_queryi(meta[j][0]->index, tid, locus, locus + 1)) == nullptr) {
            exit(EXIT_FAILURE);
        }
    }

    int max_depth = INT_MAX;
    for (int j = 0; j < num_samples; ++j) {
        bam_handler.emplace_back(bam_mplp_init(1, mplp_func_plain, (void **) (meta[j])));
        bam_mplp_init_overlaps(bam_handler.back());
        bam_mplp_set_maxcnt(bam_handler.back(), max_depth);
    }
    // begin pileup
    ReadColumns read_col;
    read_col.reserve(num_samples);
    int *tid_arr = new int[num_samples];
    int *pos_arr = new int[num_samples];
    int *n_plp_arr = new int[num_samples];
    auto plp_arr = new const bam_pileup1_t *[num_samples];
    for (int j = 0; j < num_samples; ++j) {
        std::vector<Read> reads;
        int ret;
        while ((ret = bam_mplp_auto(bam_handler[j], &tid_arr[j], &pos_arr[j], &n_plp_arr[j], &plp_arr[j])) > 0) {
            if (pos_arr[j] == locus) {

                reads.reserve(n_plp_arr[j]);
                const bam_pileup1_t *p = plp_arr[j];
                for (int i = 0; i < n_plp_arr[j]; ++i, ++p) {
                    reads.emplace_back(Read{static_cast<uint8_t >(seq_nt16_str[bam_seqi(bam_get_seq(p->b), p->qpos)]),
                                            bam_get_qual(p->b)[2]});
                }
//                std::cout << meta[j][0]->header->target_name[tid_arr[j]] << '\t' << pos_arr[j] + 1 << '\t'
//                          << n_plp_arr[j] << std::endl;
//                print_pileup_seq(plp_arr[j], n_plp_arr[j]);
//                std::cout << std::endl;
                break;
            }
        }
        read_col.emplace_back(reads);
    }
    delete[](tid_arr);
    delete[](pos_arr);
    delete[](n_plp_arr);
    delete[](plp_arr);
    return read_col;
}

std::string uint2str(const uint8_t *seq, int len) {
    std::string temp;
    for (int i = 0; i < len; ++i) {
        temp.push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }
    return temp;
}

int BamStreamer::print_pileup_seq(const bam_pileup1_t *p, int n) {
    for (int i = 0; i < n; ++i, ++p) {
        uint8_t *seq = bam_get_seq(p->b);
//            std::string temp = uint2str(seq, p->b->core.l_qseq);
        int is_rev = bam_is_rev(p->b);
        if (p->is_head) {
            std::cout << '^' << char('!' + std::min(p->b->core.qual, uint8_t(93)));
        }

        if (p->is_del) {
            std::cout << (p->is_refskip ? (is_rev ? '<' : '>') : '*');
        } else {
            const char c = seq_nt16_str[bam_seqi(seq, p->qpos)];
            std::cout << char(is_rev ? std::tolower(c) : std::toupper(c));
        }
        if (p->indel > 0) {
            int j;
            std::cout << std::setiosflags(std::ios_base::showpos) << p->indel << '(';
            for (j = 1; j <= p->indel; j++)
                std::cout << seq_nt16_str[bam_seqi(seq, p->qpos + j - p->is_del)];
            std::cout << ')';
        }
        if (p->indel < 0) {
            std::cout << std::setiosflags(std::ios_base::showpos) << p->indel << "(?)";
        }
        if (p->is_tail)
            std::cout << '$';
        std::cout << " ";
    }
    return 0;
}
