//
// Created by Chuanyi Zhang on 2019-05-21.
//

#ifndef MOSS_TYPES_H
#define MOSS_TYPES_H

#include <vector>

// internal bases are stored in 4-bit values.

namespace moss {
    enum class IUPAC_nuc : std::uint8_t {
        EQ, A, C, M, G, R, S, V, T, W, Y, H, K, D, B, N
    };

    typedef struct {
        uint8_t base;
        int qual;
    } Read;

    class Pileups {
    private:
        uint8_t ref;
        std::vector<std::vector<Read> > read_columns;
    public:
        explicit Pileups(int num_samples);

        void set_ref(char ref);

        void emplace_read_column(std::vector<Read> &col);

        const std::vector<std::vector<Read>> &get_read_columns() const;

        uint8_t get_ref() const;
    };

//    class ReadsColumn {
//    private:
//        const int n_samples;
//        std::vector<std::vector<Read>> columns;
//    public:
//        ReadsColumn(int n_reads, std::vector<Read> &column);
//
//        const int get_len() const;
//
//        const Read& get_read(int index) const;
//    };
}

#endif //MOSS_TYPES_H
