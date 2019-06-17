//
// Created by Chuanyi Zhang on 2019-05-24.
//

#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "../src/core/calling.h"

using namespace Catch;

TEST_CASE("Binomial coefficient", "[calling]") {
    REQUIRE(moss::binom(0, 0) == 1);
    REQUIRE(moss::binom(1, 0) == 1);
    REQUIRE(moss::binom(1, 1) == 1);
    REQUIRE(moss::binom(2, 0) == 1);
    REQUIRE(moss::binom(2, 1) == 2);
    REQUIRE(moss::binom(2, 2) == 1);
    REQUIRE(moss::binom(3, 0) == 1);
    REQUIRE(moss::binom(3, 1) == 3);
    REQUIRE(moss::binom(3, 2) == 3);
    REQUIRE(moss::binom(3, 3) == 1);
    REQUIRE(moss::binom(4, 0) == 1);
    REQUIRE(moss::binom(4, 1) == 4);
    REQUIRE(moss::binom(4, 2) == 6);
    REQUIRE(moss::binom(4, 3) == 4);
    REQUIRE(moss::binom(4, 4) == 1);
    // TODO: negative?
    REQUIRE(moss::binom(-1, 0) == 1);
    REQUIRE(moss::binom(10, 5) == 252);
    REQUIRE_THAT(moss::binom(100, 50), WithinAbs(1.0089134454556415e+29, 1e14));
}

TEST_CASE("Trinomial coefficient", "[calling]") {
    REQUIRE(moss::trinomial(3, 4, 3) == 4200);
    REQUIRE_THAT(moss::trinomial(33, 33, 34), WithinAbs(4.1924479242558294e+45, 1e32));
    REQUIRE(moss::trinomial(0, 0, 0) == 1);
    REQUIRE(moss::trinomial(1, 1, 1) == 6);
    REQUIRE(moss::trinomial(1, 0, 0) == 1);
    REQUIRE(moss::trinomial(0, 1, 0) == 1);
    REQUIRE(moss::trinomial(0, 0, 1) == 1);
    REQUIRE(moss::trinomial(1, 1, 0) == 2);
    REQUIRE(moss::trinomial(1, 0, 1) == 2);
    REQUIRE(moss::trinomial(0, 1, 1) == 2);
}

TEST_CASE("q-phred to probability", "[calling]") {
    REQUIRE_THAT(moss::qphred2prob(30), WithinAbs(1e-3, 1e-21));
}

TEST_CASE("Log sum exp", "[calling]") {
    std::vector<double> v1{log(0.01), log(0.02), log(0.03), log(0.04)};
    REQUIRE_THAT(moss::log_sum_exp(v1), WithinAbs(log(0.1), 1e-21));

    std::vector<double> v2{log(1e-31), log(2e-31), log(3e-31), log(4e-31)};
    REQUIRE_THAT(moss::log_sum_exp(v2), WithinAbs(log(1e-30), 1e-21));

    std::vector<double> v3{-334, -334 + log(2), -334 + log(3), -334 + log(4)};
    REQUIRE_THAT(moss::log_sum_exp(v3), WithinAbs(-334 + log(10), 1e-21));
}
