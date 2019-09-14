//
// Created by Chuanyi Zhang on 2019-05-24.
//

#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "../src/core/calling.h"

using namespace Catch;

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
    double max_outer, accu_outer;
    moss::log_sum_exp_init(max_outer, accu_outer);
    for (int i = 0; i < 2; ++i) {
        double max_item, accu;
        moss::log_sum_exp_init(max_item, accu);
        for (const auto &item : v3) {
            moss::log_sum_exp_iter(max_item, accu, item);
        }
        REQUIRE_THAT(moss::log_sum_exp_final(max_item, accu), WithinAbs(-334 + log(10), 1e-21));
        moss::log_sum_exp_iter(max_outer, accu_outer, moss::log_sum_exp_final(max_item, accu));
    }
    REQUIRE_THAT(moss::log_sum_exp_final(max_outer, accu_outer), WithinAbs(-334 + log(20), 1e-21));
}
