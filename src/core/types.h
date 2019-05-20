//
// Created by Chuanyi Zhang on 2019-05-21.
//

#ifndef MOSS_TYPES_H
#define MOSS_TYPES_H

#include <vector>

namespace moss {
    typedef struct {
        uint8_t base;
        int qual;
    } Read;

    using ReadColumns = std::vector<std::vector<Read>>;

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
