//
// Created by Chuanyi Zhang on 2019-05-16.
//

#include "bam_io.h"
#include <iostream>
#include <climits>
#include <iomanip>

using namespace moss;

BamStreamer::BamStreamer(std::string ref_file_name, std::vector<std::string> bam_file_names, MapContigLoci loci,
                         Stepper stepper, int min_baseQ) : num_samples(bam_file_names.size()), min_base_qual(min_baseQ),
                                                           reference(std::move(ref_file_name)), loci(loci),
                                                           tids(bam_file_names.size(), -1) {
    switch (stepper) {
        case Stepper::nofilter:
            pileup_func = BamStreamer::mplp_func_plain;
            break;
        case Stepper::samtools:
        default:
            pileup_func = BamStreamer::mplp_func_samtools;
            break;
    }
    ref_fp = fai_load(reference.c_str());
    meta.reserve(num_samples);
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
    regions.reserve(num_samples);
    windows.reserve(num_samples);
    iters.reserve(num_samples);
    actives.reserve(num_samples);
    buffers.reserve(num_samples);
    for (int idx_sample = 0; idx_sample < num_samples; ++idx_sample) {
        regions.emplace_back(std::vector<hts_reglist_t>());
        regions.back().reserve(loci.size());
        windows.emplace_back(std::make_pair(0, 0));
        actives.emplace_back(std::deque<locus_t>());
        buffers.emplace_back(Buffer());
    }
    for (const auto &contig : this->loci) {
        auto intervals = new hts_pair32_t[contig.second.size()];
        // TODO: avoid copying this?
        char *reg = new char[contig.first.size() + 1];
        std::strcpy(reg, contig.first.c_str());
        int idx_loc = 0;
        uint32_t min = 0xFFFFFFFF,
            max = 0;
        for (const auto &locus : contig.second) {
            intervals[idx_loc].beg = locus;
            intervals[idx_loc].end = locus + 1;
            min = locus < min ? locus : min;
            max = locus + 1 > max ? locus + 1 : max;
            idx_loc++;
        }
        for (int idx_sample = 0; idx_sample < num_samples; ++idx_sample) {
            hts_reglist_t list{
                reg,
                bam_name2id(meta[idx_sample][0]->header, contig.first.c_str()),
                intervals,
                static_cast<uint32_t>(contig.second.size()),
                min,
                max
            };
            regions[idx_sample].emplace_back(list);
        }
    }
    for (int j = 0; j < num_samples; ++j) {
        if ((meta[j][0]->iter = sam_itr_regions(meta[j][0]->index, meta[j][0]->header, regions[j].data(),
                                                regions[j].size())) == nullptr) {
            exit(EXIT_FAILURE);
        }
    }
}

BamStreamer::~BamStreamer() {
    fai_destroy(ref_fp);
    for (auto &item : meta) {
        sam_close(item[0]->sam_fp);
        bam_hdr_destroy(item[0]->header);
        delete[](item[0]);
        delete[](item);
    }
    for (const auto &reg : regions[0]) {
        delete[](reg.intervals);
        delete[](reg.reg);
    }
}

int BamStreamer::mplp_func_plain(void *data, bam1_t *b) {
    int ret;
    auto d = static_cast<data_t *>(data);
    return ret;
}

int BamStreamer::mplp_func_samtools(void *data, bam1_t *b) {
    int ret,
            skip = 0;
    auto d = static_cast<data_t *>(data);
    do {
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

// TODO: do this contig-wise, make sure all samples are in the same contig simultaneously
Pileups BamStreamer::get_column() {
    Pileups read_col(num_samples);

    // begin pileup
    bam1_t *b = bam_init1();
    int ret;
    /*!
     * \details For each sample, get reads until the first active locus in the queue is finished,
     * then emplace in the output vector
     */
    for (int j = 0; j < num_samples; ++j) {
        int cnt = 0;
        while ((ret = sam_itr_multi_next(meta[j][0]->sam_fp, meta[j][0]->iter, b)) >= 0) {
            uint32_t *cigar = bam_get_cigar(b);
//            for (int i = 0; i < b->core.n_cigar; ++i) {
//                std::cout << (cigar[i] >> BAM_CIGAR_SHIFT) << BAM_CIGAR_STR[cigar[i] & BAM_CIGAR_MASK];
//            }
            locus_t begin = b->core.pos;
            locus_t end = b->core.pos + b->core.l_qseq - 1;
//            std::cout << '\t' << b->core.pos << '\t' << b->core.l_qseq << '\t' << end << std::endl;
            cnt++;
            if (b->core.tid < 0 || (b->core.flag & BAM_FUNMAP)) {
                continue;
            }
            // update window
            windows[j].first = b->core.pos;
            if (end > windows[j].second) {
                windows[j].second = end;
            }
            // check iter
            if (b->core.tid != tids[j]) {
                tids[j] = b->core.tid;
                if (iters.size() <= j) {
                    iters.push_back(this->loci.at(std::string(meta[j][0]->header->target_name[b->core.tid])).cbegin());
                } else {
                    iters[j] = (this->loci.at(std::string(meta[j][0]->header->target_name[b->core.tid])).cbegin());
                }
            }
            // find new active loci
            while (*iters[j] <= windows[j].second) {
                actives[j].emplace_back(*iters[j]);
                buffers[j].emplace_back(std::vector<Read>());
                iters[j]++;
            }
            // TODO: resolve cigar
            // push into queue
            int idx_pos = 0;
            for (const auto &pos : actives[j]) {
                int qpos = pos - windows[j].first;          // use cigar!
                if (qpos >= 0)
                {
                    uint8_t base_qual = bam_get_qual(b)[qpos];
                    if (base_qual >= min_base_qual) {
                        buffers[j][idx_pos].push_back(Read{static_cast<uint8_t >(bam_seqi(bam_get_seq(b), qpos)),
                                                           base_qual});
                    }
                }
                idx_pos++;
            }
            // Until first in the queue finished
            if (windows[j].first > actives[j][0]) {
                break;
            }
        }
        // find reference, only once
        if (j == 0) {
            int len_seq;
            char *temp = faidx_fetch_seq(ref_fp, meta[j][0]->header->target_name[b->core.tid], actives[j].front(),
                                         actives[j].front(), &len_seq);
            read_col.set_ref(temp[0]);
            free(temp);
        }

        actives[j].pop_front();
        read_col.emplace_read_column(std::move(buffers[j].front()));
        buffers[j].pop_front();
//        std::cout << cnt << std::endl;
    }

    bam_destroy1(b);
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

//typedef struct {
//    int k, x, y, end;
//} cstate_t;
//
///* s->k: the index of the CIGAR operator that has just been processed.
//   s->x: the reference coordinate of the start of s->k
//   s->y: the query coordiante of the start of s->k
// */
//int resolve_cigar2(bam1_t *b, int32_t pos, cstate_t *s)
//{
//#define _cop(c) ((c)&BAM_CIGAR_MASK)
//#define _cln(c) ((c)>>BAM_CIGAR_SHIFT)
//
////    bam1_t *b = p->b;
//    bam1_core_t *c = &b->core;
//    uint32_t *cigar = bam_get_cigar(b);
//    int k;
//    // determine the current CIGAR operation
//    //fprintf(stderr, "%s\tpos=%d\tend=%d\t(%d,%d,%d)\n", bam_get_qname(b), pos, s->end, s->k, s->x, s->y);
//    if (s->k == -1) { // never processed
//        p->qpos = 0;
//        if (c->n_cigar == 1) { // just one operation, save a loop
//            if (_cop(cigar[0]) == BAM_CMATCH || _cop(cigar[0]) == BAM_CEQUAL || _cop(cigar[0]) == BAM_CDIFF) s->k = 0, s->x = c->pos, s->y = 0;
//        } else { // find the first match or deletion
//            for (k = 0, s->x = c->pos, s->y = 0; k < c->n_cigar; ++k) {
//                int op = _cop(cigar[k]);
//                int l = _cln(cigar[k]);
//                if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
//                    op == BAM_CEQUAL || op == BAM_CDIFF) break;
//                else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
//            }
//            assert(k < c->n_cigar);
//            s->k = k;
//        }
//    } else { // the read has been processed before
//        int op, l = _cln(cigar[s->k]);
//        if (pos - s->x >= l) { // jump to the next operation
//            assert(s->k < c->n_cigar); // otherwise a bug: this function should not be called in this case
//            op = _cop(cigar[s->k+1]);
//            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) { // jump to the next without a loop
//                if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
//                s->x += l;
//                ++s->k;
//            } else { // find the next M/D/N/=/X
//                if (_cop(cigar[s->k]) == BAM_CMATCH|| _cop(cigar[s->k]) == BAM_CEQUAL || _cop(cigar[s->k]) == BAM_CDIFF) s->y += l;
//                s->x += l;
//                for (k = s->k + 1; k < c->n_cigar; ++k) {
//                    op = _cop(cigar[k]), l = _cln(cigar[k]);
//                    if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP || op == BAM_CEQUAL || op == BAM_CDIFF) break;
//                    else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) s->y += l;
//                }
//                s->k = k;
//            }
//            assert(s->k < c->n_cigar); // otherwise a bug
//        } // else, do nothing
//    }
//    { // collect pileup information
//        int op, l;
//        op = _cop(cigar[s->k]); l = _cln(cigar[s->k]);
//        p->is_del = p->indel = p->is_refskip = 0;
//        if (s->x + l - 1 == pos && s->k + 1 < c->n_cigar) { // peek the next operation
//            int op2 = _cop(cigar[s->k+1]);
//            int l2 = _cln(cigar[s->k+1]);
//            if (op2 == BAM_CDEL) p->indel = -(int)l2;
//            else if (op2 == BAM_CINS) p->indel = l2;
//            else if (op2 == BAM_CPAD && s->k + 2 < c->n_cigar) { // no working for adjacent padding
//                int l3 = 0;
//                for (k = s->k + 2; k < c->n_cigar; ++k) {
//                    op2 = _cop(cigar[k]); l2 = _cln(cigar[k]);
//                    if (op2 == BAM_CINS) l3 += l2;
//                    else if (op2 == BAM_CDEL || op2 == BAM_CMATCH || op2 == BAM_CREF_SKIP || op2 == BAM_CEQUAL || op2 == BAM_CDIFF) break;
//                }
//                if (l3 > 0) p->indel = l3;
//            }
//        }
//        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
//            p->qpos = s->y + (pos - s->x);
//        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
//            p->is_del = 1; p->qpos = s->y; // FIXME: distinguish D and N!!!!!
//            p->is_refskip = (op == BAM_CREF_SKIP);
//        } // cannot be other operations; otherwise a bug
//        p->is_head = (pos == c->pos); p->is_tail = (pos == s->end);
//    }
//    p->cigar_ind = s->k;
//    return 1;
//}
