//
// Created by Chuanyi Zhang on 2019-05-24.
//

#include "catch.hpp"
#include "../src/core/types.h"
#include "htslib/hts.h"

using namespace Catch;

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

TEST_CASE("BaseSet", "[types]") {
    using moss::operator"" _8;
    auto set = moss::BaseSet(static_cast<uint8_t>(moss::IUPAC_nuc::M));

    REQUIRE(set.size() == 2);
    REQUIRE(set.contain(static_cast<uint8_t>(moss::IUPAC_nuc::A)));
    REQUIRE(set.contain(static_cast<uint8_t>(moss::IUPAC_nuc::C)));

    REQUIRE(moss::BaseSet(0x00_8).size() == 0);
    REQUIRE(moss::BaseSet(0x01_8).size() == 1);
    REQUIRE(moss::BaseSet(0x02_8).size() == 1);
    REQUIRE(moss::BaseSet(0x03_8).size() == 2);
    REQUIRE(moss::BaseSet(0x04_8).size() == 1);
    REQUIRE(moss::BaseSet(0x05_8).size() == 2);
    REQUIRE(moss::BaseSet(0x06_8).size() == 2);
    REQUIRE(moss::BaseSet(0x07_8).size() == 3);
    REQUIRE(moss::BaseSet(0x08_8).size() == 1);
    REQUIRE(moss::BaseSet(0x09_8).size() == 2);
    REQUIRE(moss::BaseSet(0x0A_8).size() == 2);
    REQUIRE(moss::BaseSet(0x0B_8).size() == 3);
    REQUIRE(moss::BaseSet(0x0C_8).size() == 2);
    REQUIRE(moss::BaseSet(0x0D_8).size() == 3);
    REQUIRE(moss::BaseSet(0x0E_8).size() == 3);
    REQUIRE(moss::BaseSet(0x0F_8).size() == 4);

    auto difference = moss::BaseSet::set_difference(0x0c_8, 0x0a_8);
    REQUIRE(difference.get_set() == 0x04_8);
    REQUIRE_THAT(difference.get_base_list(), Equals(std::vector<uint8_t>({0x04_8})));

    auto comple = moss::BaseSet(0x06_8).complement();
    REQUIRE(comple.get_set() == 0x09_8);
    REQUIRE_THAT(comple.get_base_list(), Equals(std::vector<uint8_t>({0x01_8, 0x08_8})));
}
