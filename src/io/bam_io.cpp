//
// Created by Chuanyi Zhang on 2019-05-16.
//

#include "bam_io.h"
#include <iostream>
#include <climits>
#include <iomanip>

using namespace moss;

BamStreamer::BamStreamer(std::string ref_file_name, std::vector<std::string> bam_file_names, Stepper stepper,
                         int min_baseQ) : reference(std::move(ref_file_name)), min_base_qual(min_baseQ) {
    switch (stepper) {
        case Stepper::nofilter:
            pileup_func = BamStreamer::mplp_func_plain;
            break;
        case Stepper::samtools:
        default:
            pileup_func = BamStreamer::mplp_func_samtools;
            break;
    }
    min_base_qual = min_baseQ;
    num_samples = bam_file_names.size();
    ref_fp = fai_load(reference.c_str());
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
    fai_destroy(ref_fp);
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

int BamStreamer::mplp_func_plain(void *data, bam1_t *b) {
    int ret;
    auto d = static_cast<data_t *>(data);
    ret = sam_itr_next(d->sam_fp, d->iter, b);
    return ret;
}

int BamStreamer::mplp_func_samtools(void *data, bam1_t *b) {
    int ret,
            skip = 0;
    auto d = static_cast<data_t *>(data);
    do {
        ret = sam_itr_next(d->sam_fp, d->iter, b);
        if (ret < 0) {
            break;
        }
        if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP)) {
            skip = 1;
            continue;
        }

//        if (ma->conf->rflag_require && !(ma->conf->rflag_require&b->core.flag)) { skip = 1; continue; }
//        if (ma->conf->rflag_filter && ma->conf->rflag_filter&b->core.flag) { skip = 1; continue; }
//        if (ma->conf->bed && ma->conf->all == 0) { // test overlap
//            skip = !bed_overlap(ma->conf->bed, ma->h->target_name[b->core.tid], b->core.pos, bam_endpos(b));
//            if (skip) continue;
//        }
//        if (ma->conf->rghash) { // exclude read groups
//            uint8_t *rg = bam_aux_get(b, "RG");
//            skip = (rg && khash_str2int_get(ma->conf->rghash, (const char*)(rg+1), NULL)==0);
//            if (skip) continue;
//        }
//        if (ma->conf->flag & MPLP_ILLUMINA13) {
//            int i;
//            uint8_t *qual = bam_get_qual(b);
//            for (i = 0; i < b->core.l_qseq; ++i)
//                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
//        }
//
//        if (ma->conf->fai && b->core.tid >= 0) {
//            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
//            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
//                fprintf(stderr,"[%s] Skipping because %d is outside of %d [ref:%d]\n",
//                        __func__, b->core.pos, ref_len, b->core.tid);
//                skip = 1;
//                continue;
//            }
//        } else {
//            has_ref = 0;
//        }
//
//        skip = 0;
//        if (has_ref && (ma->conf->flag&MPLP_REALN)) sam_prob_realn(b, ref, ref_len, (ma->conf->flag & MPLP_REDO_BAQ)? 7 : 3);
//        if (has_ref && ma->conf->capQ_thres > 10) {
//            int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
//            if (q < 0) skip = 1;
//            else if (b->core.qual > q) b->core.qual = q;
//        }
//        if (b->core.qual < ma->conf->min_mq) skip = 1;
//        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) skip = 1;
    } while (skip);
    return ret;
}

Pileups BamStreamer::get_column(std::string contig, int locus) {
    Pileups read_col(num_samples);
    for (int j = 0; j < num_samples; ++j) {
        int tid = bam_name2id(meta[j][0]->header, contig.c_str());
        if ((meta[j][0]->iter = sam_itr_queryi(meta[j][0]->index, tid, locus, locus + 1)) == nullptr) {
            exit(EXIT_FAILURE);
        }
    }

    int max_depth = INT_MAX;
    for (int j = 0; j < num_samples; ++j) {
        bam_handler.emplace_back(bam_mplp_init(1, pileup_func, (void **) (meta[j])));
        bam_mplp_init_overlaps(bam_handler.back());
        bam_mplp_set_maxcnt(bam_handler.back(), max_depth);
    }
    // find reference
    int len_seq;
    char *temp = faidx_fetch_seq(ref_fp, contig.c_str(), locus, locus, &len_seq);
    read_col.set_ref(temp[0]);
    free(temp);
    // begin pileup
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
                    uint8_t base_qual = bam_get_qual(p->b)[p->qpos];
                    if (base_qual >= min_base_qual) {
                        reads.emplace_back(
                                Read{static_cast<uint8_t >(bam_seqi(bam_get_seq(p->b), p->qpos)),
                                     base_qual});
                    }
                }
                break;
            }
        }
        read_col.emplace_read_column(reads);
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
