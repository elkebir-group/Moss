//
// Created by Chuanyi Zhang on 2019-05-24.
//

#include <htslib/hts.h>
#include "types.h"

using namespace moss;

Pileups::Pileups(int num_samples) {
    read_columns.reserve(num_samples);
}

void Pileups::set_ref(char ref) {
    Pileups::ref = seq_nt16_table[ref];
}

// TODO: rvalue ref here?
void Pileups::emplace_read_column(std::vector<Read> &&col) {
    read_columns.emplace_back(col);
}

const std::vector<std::vector<Read>> &Pileups::get_read_columns() const {
    return read_columns;
}

uint8_t Pileups::get_ref() const {
    return ref;
}

BaseSet::BaseSet(uint8_t base) : set(base) {
    base_list.reserve(count_1bits[base]);
    for (int i = 0; i < 4; ++i) {
        if ((base & (1 << i)) != 0) {
            base_list.push_back(static_cast<uint8_t >(1 << i));
        }
    }
}

uint8_t BaseSet::get_set(){
    return set;
}

const unsigned BaseSet::count_1bits[16] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

const std::vector<uint8_t> &BaseSet::get_base_list() {
    return base_list;
}

void BaseSet::add_base(uint8_t base) {
//    set |= seq_nt16_table[base];
    set |= base;
}
