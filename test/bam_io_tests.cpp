//
// Created by Chuanyi Zhang on 2019-06-18.
//

#include "catch.hpp"
#include "../src/io/bam_io.cpp"

using namespace Catch;

TEST_CASE("CIGAR resolving", "[bam]") {
    bam1_core_t core{};
    core.n_cigar = 3;
    core.pos = 15;
    core.l_qname = 8;
    core.l_qseq = 11;
    core.l_extranul = 4;

    uint8_t data[37]{};
    uint32_t cigar[3]{
        (6 << BAM_CIGAR_SHIFT) | BAM_CMATCH,
        (14 << BAM_CIGAR_SHIFT) | BAM_CREF_SKIP,
        (5 << BAM_CIGAR_SHIFT) | BAM_CMATCH
    };
    uint8_t seq[6]{(1 << 4) | 8, (1 << 4) | 4, (2 << 4) | 8, (8 << 4) | 2, (1 << 4) | 4, (2 << 4)}; // ATAGCTTCAGC

#if HTS_VERSION >= 101000
    bam1_t b{
        core,       // core
        0,          // id
        data,
        static_cast<int>(core.l_qname + core.n_cigar * 4 + (core.l_qseq + 1) / 2 + core.l_qseq),         // l_data
        2048,       // m_data
        0,
    };
#else
    bam1_t b{
        core,       // core
        static_cast<int>(core.l_qname + core.n_cigar * 4 + (core.l_qseq + 1) / 2 + core.l_qseq),         // l_data
        2048,       // m_data
        data,
    };
#endif

    std::strcpy(bam_get_qname(&b), "r004");
    std::memcpy(bam_get_cigar(&b), cigar, b.core.n_cigar * 4);
    std::memcpy(bam_get_seq(&b), seq, sizeof(seq));
    std::memset(bam_get_qual(&b), 40, b.core.l_qseq);

    PileupMeta p{};
    resolve_cigar(&b, 39, &p);
    REQUIRE(p.qpos == 10);
    resolve_cigar(&b, 20, &p);
    REQUIRE(p.qpos == 5);
    resolve_cigar(&b, 24, &p);
    REQUIRE(p.is_refskip == 1);
    resolve_cigar(&b, 35, &p);
    REQUIRE(p.qpos == 6);
    resolve_cigar(&b, 39, &p);
    REQUIRE(p.qpos == 10);
}
