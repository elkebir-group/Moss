//
// Created by Chuanyi Zhang on 2019-05-24.
//

#include <htslib/hts.h>
#include "types.h"

moss::Pileups::Pileups(int num_samples) {
    read_columns.reserve(num_samples);
}

void moss::Pileups::set_ref(char ref) {
    Pileups::ref = seq_nt16_table[ref];
}

void moss::Pileups::emplace_read_column(std::vector<moss::Read> &col) {
    read_columns.emplace_back(col);
}

const std::vector<std::vector<moss::Read>> &moss::Pileups::get_read_columns() const {
    return read_columns;
}

uint8_t moss::Pileups::get_ref() const {
    return ref;
}
