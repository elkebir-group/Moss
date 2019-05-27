//
// Created by Chuanyi Zhang on 2019-05-24.
//

#include "catch.hpp"
#include "../src/core/types.h"
#include "htslib/hts.h"

TEST_CASE("IUPAC notation defined", "[types]") {
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::EQ) == 0);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::A) == 1);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::C) == 2);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::M) == 3);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::G) == 4);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::R) == 5);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::S) == 6);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::V) == 7);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::T) == 8);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::W) == 9);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::Y) == 10);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::H) == 11);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::K) == 12);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::D) == 13);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::B) == 14);
    REQUIRE(static_cast<uint8_t>(moss::IUPAC_nuc::N) == 15);
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::EQ)] == '=');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::A)] == 'A');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::C)] == 'C');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::M)] == 'M');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::G)] == 'G');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::R)] == 'R');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::S)] == 'S');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::V)] == 'V');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::T)] == 'T');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::W)] == 'W');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::Y)] == 'Y');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::H)] == 'H');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::K)] == 'K');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::D)] == 'D');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::B)] == 'B');
    REQUIRE(seq_nt16_str[static_cast<uint8_t>(moss::IUPAC_nuc::N)] == 'N');
}
