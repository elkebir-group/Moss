//
// Created by Chuanyi Zhang on 2019-05-16.
//

#include "bam_io.h"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>

using namespace moss;

typedef struct {
    int32_t qpos;
    int indel, level;
    uint32_t is_del: 1, is_head: 1, is_tail: 1, is_refskip: 1, aux: 28;
    int cigar_ind;
} PileupMeta;

void resolve_cigar(bam1_t *b, int64_t pos, PileupMeta *p) {
    bam1_core_t *c = &b->core;
    uint32_t *cigar = bam_get_cigar(b);
    int64_t dist = pos - c->pos;
    int qpos = 0;
    bool is_finish = false;
    for (int k = 0; k < c->n_cigar; k++) {
        int op = bam_cigar_op(cigar[k]);
        int l = bam_cigar_oplen(cigar[k]);
        switch (op) {
            case BAM_CMATCH:
            case BAM_CEQUAL:
            case BAM_CDIFF:
                if (dist <= l - 1) {
                    p->qpos = qpos + dist;
                    is_finish = true;
                    break;
                }
                dist -= l;
                qpos += l;
                break;
            case BAM_CINS:
            case BAM_CSOFT_CLIP:
                qpos += l;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                if (dist <= l - 1) {
                    p->is_del = 1;
                    p->is_refskip = (op == BAM_CREF_SKIP);
                    is_finish = true;
                    break;
                }
                dist -= l;
                break;
            case BAM_CHARD_CLIP:
            case BAM_CPAD:
                break;
        }
        if (is_finish) {
            break;
        }
    }
}

SingleBamStreamer::SingleBamStreamer(const std::string &ref_file_name,
                                     const std::string &bam_file_name,
                                     MapContigLoci loci,
                                     int min_baseQ,
                                     int min_mapQ,
                                     bool filter_edit_distance)
    : loci(std::move(loci)),
      min_base_qual(min_baseQ),
      min_map_qual(min_mapQ),
      tid(-1),
      is_filter_edit_distance(filter_edit_distance) {
    htsFormat fmt = {};
    samFile *fp = sam_open_format(bam_file_name.c_str(), "rb", &fmt);
    if (fp == nullptr) {
        throw std::runtime_error("Failed to open BAM file " + bam_file_name);
    }
    hts_idx_t *idx = sam_index_load(fp, bam_file_name.c_str());
    if (idx == nullptr) {
        throw std::runtime_error("Failed to open BAM index of " + bam_file_name);
    }
    bam_hdr_t *hdr = sam_hdr_read(fp);
    if (hdr == nullptr) {
        throw std::runtime_error(
            "Failed to fetch header from BAM file " + bam_file_name + ".fai , the file may be broken.");
    }
    auto tmp = new data_t *[1];
    tmp[0] = new data_t[1];
    tmp[0][0] = {fp, nullptr, idx, hdr};
    meta = tmp;
    region = std::vector<hts_reglist_t>();
    region.reserve(loci.size());
    window = std::make_pair(0, 0);
    actives = std::deque<locus_t>();
    buffer = Buffer();
    for (const auto &contig : this->loci) {
        auto intervals = new hts_pair32_t[contig.second.size()];
        // TODO: avoid copying this?
        char *reg = new char[contig.first.size() + 1];
        std::strcpy(reg, contig.first.c_str());
        int idx_loc = 0;
        uint32_t min = 0xFFFFFFFF,
            max = 0;
        for (const auto &l : contig.second) {
            intervals[idx_loc].beg = l.first;
            intervals[idx_loc].end = l.first + 1;
            min = l.first < min ? l.first : min;
            max = l.first + 1 > max ? l.first + 1 : max;
            idx_loc++;
        }
#if HTS_VERSION >= 101000
        hts_reglist_t list{
            reg,
            intervals,
            bam_name2id(meta[0]->header, contig.first.c_str()),
            static_cast<uint32_t>(contig.second.size()),
            min,
            max
        };
#else
        hts_reglist_t list{
            reg,
            bam_name2id(meta[0]->header, contig.first.c_str()),
            intervals,
            static_cast<uint32_t>(contig.second.size()),
            min,
            max
        };
#endif
        region.emplace_back(list);
    }
    if ((meta[0]->iter = sam_itr_regions(meta[0]->index, meta[0]->header, region.data(),
                                         region.size())) == nullptr) {
        exit(EXIT_FAILURE);
    }
}

SingleBamStreamer::~SingleBamStreamer() {
    sam_close(meta[0]->sam_fp);
    bam_hdr_destroy(meta[0]->header);
    delete[](meta[0]);
    delete[](meta);
    delete[](region[0].intervals);
    delete[](region[0].reg);
}

std::vector<Read> SingleBamStreamer::get_column(std::string &ret_contig, int &ret_pos) {
    bool not_found = true;
    /// begin pileup
    bam1_t *read = bam_init1();
    int ret;
    /*!
     * \details For each sample, get reads until the first active locus in the queue is finished,
     * then emplace in the output vector
     */
    int cnt = 0;
    if (actives.empty() || window.first <= actives[0]) {
        do {
            ret = sam_itr_multi_next(meta[0]->sam_fp, meta[0]->iter, read);
            if (ret < 0) {
                break;
            }
            locus_t begin = read->core.pos;
            locus_t end = bam_endpos(read);
            cnt++;
            // next contig
            if (read->core.tid != tid) {
                tid = read->core.tid;
                window.second = 0;
                iter = this->loci.at(std::string(meta[0]->header->target_name[read->core.tid])).cbegin();
            }
            if (read->core.tid < 0 || (read->core.flag & FAIL_FLAGS)) {
                continue;
            }
            // TODO: optional
            if (read->core.qual < min_map_qual) {
                continue;
            }
            if (is_filter_edit_distance) {
                uint8_t *nm = bam_aux_get(read, "NM");
                if (nm) {
                    int64_t edit_dist = bam_aux2i(nm);
                    if (edit_dist > 3) {
                        continue;
                    }
                }
            }
            /// update window
            window.first = begin;
            if (end > window.second) {
                window.second = end;
            }

            /// find new active loci
            while (iter->first <= window.second) {
                if (iter !=
                    this->loci.at(std::string(meta[0]->header->target_name[read->core.tid])).cend()) {
                    actives.emplace_back(iter->first);
                    buffer.emplace_back(std::vector<Read>());
                    ++iter;
                } else {
                    break;
                }
            }
            /// push into queue
            int idx_pos = 0;
            for (const auto &pos : actives) {
                if (pos >= window.first && pos < end) {
                    PileupMeta p{};
                    resolve_cigar(read, pos, &p);
                    if (!p.is_del && p.qpos >= 0 && p.qpos < read->core.l_qseq) {
                        uint8_t base_qual = bam_get_qual(read)[p.qpos];
                        if (base_qual >= min_base_qual) {
                            buffer[idx_pos].emplace_back(static_cast<uint8_t>(bam_seqi(bam_get_seq(read), p.qpos)),
                                                         base_qual,
                                                         static_cast<bool>(read->core.flag & BAM_FREVERSE),
                                                         std::string(bam_get_qname(read)));
                        }
                    }
                }
                idx_pos++;
            }
        } while (window.first <= actives[0]); // FIXME: Conditional jump or move depends on uninitialised value(s)
    }

    if (!actives.empty()) {
        ret_contig = meta[0]->header->target_name[tid];
        ret_pos = actives.front();
        actives.pop_front();
    }

    std::vector<Read> sample_column;
    if (buffer.empty()) {
        sample_column = std::vector<Read>();
    } else {
        sample_column = std::move(buffer.front());
        buffer.pop_front();
    }

    bam_destroy1(read);
    return sample_column;
}

PairedBamStreamer::PairedBamStreamer(const std::string &ref_file_name,
                                     const std::string &original_bam_file_name,
                                     const std::string &realigned_bam_file_name,
                                     const MapContigLoci loci,
                                     int min_baseQ,
                                     int min_mapQ,
                                     bool filter_edit_distance)
    : original(ref_file_name, original_bam_file_name, loci, min_baseQ, min_mapQ, filter_edit_distance),
      realigned(ref_file_name, realigned_bam_file_name, loci, min_baseQ, min_mapQ, filter_edit_distance) {}

std::vector<Read> PairedBamStreamer::get_column(std::string &ret_contig, int &ret_pos) {
    std::string contig;
    int pos;
    std::vector<Read> orig_col = original.get_column(ret_contig, ret_pos);
    std::vector<Read> realign_col = realigned.get_column(contig, pos);
    if (realign_col.empty()) {
        return orig_col;
    } else {
        for (const auto &realigned_read : realign_col) {
            bool is_new = true;
            for (auto &&original_read : orig_col) {
                if (original_read.name == realigned_read.name) {
                    original_read.base = realigned_read.base;
                    original_read.qual = realigned_read.qual;
                    is_new = false;
                    break;
                }
            }
            if (is_new) {
                orig_col.push_back(realigned_read);
            }
        }
    }
    return orig_col;
}

MultiBamStreamer::MultiBamStreamer(std::string ref_file_name, const std::vector<std::string> &bam_file_names,
                                   const std::vector<std::string> &realigned_file_names,
                                   const MapContigLoci &loci, int min_baseQ, int min_mapQ, bool filter_edit_distance)
    : num_samples(bam_file_names.size()),
      min_base_qual(min_baseQ),
      min_map_qual(min_mapQ),
      reference(std::move(ref_file_name)),
      loci(loci),
      is_filter_edit_distance(filter_edit_distance) {
    // ref_fp = fai_load(reference.c_str());
    ref_fp = fai_load3(reference.c_str(), (reference + ".fai").c_str(), nullptr, 0x0);
    if (ref_fp == nullptr) {
        throw std::runtime_error("Failed to open reference file " + reference);
    }
    if (bam_file_names.size() != realigned_file_names.size()) {
        throw std::runtime_error("Different numbers of original and realigned BAMs in vectors");
    }
    streams.reserve(num_samples);
    for (int idx = 0; idx < bam_file_names.size(); idx++) {
        if (realigned_file_names[idx].empty()) {
            streams.push_back(std::move(std::unique_ptr<SingleBamStreamer>(
                new SingleBamStreamer(ref_file_name,
                                      bam_file_names[idx],
                                      loci,
                                      min_baseQ,
                                      min_mapQ,
                                      filter_edit_distance))));
        } else {
            streams.push_back(std::move(std::unique_ptr<PairedBamStreamer>(
                new PairedBamStreamer(ref_file_name,
                                      bam_file_names[idx],
                                      realigned_file_names[idx],
                                      loci,
                                      min_baseQ,
                                      min_mapQ,
                                      filter_edit_distance))));
        }
    }
}

MultiBamStreamer::~MultiBamStreamer() {
    fai_destroy(ref_fp);
}

// TODO: do this contig-wise, make sure all samples are in the same contig simultaneously
Pileups MultiBamStreamer::get_column() {
    Pileups read_col(num_samples);
    bool not_found = true;
    std::string contig;
    int pos;
    /*!
     * \details For each sample, get reads until the first active locus in the queue is finished,
     * then emplace in the output vector
     */
    for (int j = 0; j < num_samples; ++j) {
        read_col.emplace_read_column(streams[j]->get_column(contig, pos));

        if (not_found && contig.length() > 0) {
            int len_seq;
            char *temp = faidx_fetch_seq(ref_fp, contig.c_str(), pos, pos, &len_seq);
            if (temp != nullptr) {
                read_col.set_ref(temp[0]);
                free(temp);
            } else {
                throw std::runtime_error("MultiBamStreamer::get_column: reference error.");
            }
            not_found = false;
        }
    }
    return read_col;
}

std::string uint2str(const uint8_t *seq, int len) {
    std::string temp;
    for (int i = 0; i < len; ++i) {
        temp.push_back(seq_nt16_str[bam_seqi(seq, i)]);
    }
    return temp;
}
