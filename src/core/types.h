//
// Created by Chuanyi Zhang on 2019-05-21.
//

#ifndef MOSS_TYPES_H
#define MOSS_TYPES_H

#include <vector>
#include <set>
#include <map>

// internal bases are stored in 4-bit values.

namespace moss {
    using locus_t = unsigned long;

    // TODO: is std::set really needed here?
    using MapContigLoci = std::map<std::string, std::set<locus_t>>;

    enum class IUPAC_nuc : std::uint8_t {
        EQ, A, C, M, G, R, S, V, T, W, Y, H, K, D, B, N
    };

    inline std::uint8_t operator "" _8(unsigned long long value) {
        return static_cast<std::uint8_t>(value);
    }

    class BaseSet {
    private:
        uint8_t set;
        static const unsigned count_1bits[16];
    public:
        BaseSet() = default;

        explicit BaseSet(uint8_t base);

        unsigned size();

        bool contain(uint8_t nuc);

        bool is_valid();

        uint8_t get_set();

        static BaseSet set_difference(uint8_t a, uint8_t b);

        BaseSet complement();

        void add_base(uint8_t base);

        class iter {
        public:
            iter(uint8_t pos, uint8_t set) : pos(pos), set(set) {}

            iter operator++() {
                do {
                    pos <<= 1;
                } while (((pos & set) == 0) && (pos < (1 << 4)));
                return *this;
            }

            bool operator!=(const iter &other) const { return pos != other.pos; }

            const uint8_t operator*() const { return pos; }

        private:
            uint8_t pos;
            uint8_t set;
        };

        iter begin() const {
            uint8_t mask = 1;
            while ((mask & set) == 0) { mask <<= 1; }
            return {mask, set};
        }

        iter end() const { return {16, set}; }
    };

    inline unsigned moss::BaseSet::size() {
        return count_1bits[set];
    }

    inline bool BaseSet::contain(uint8_t nuc) {
        return (set & nuc) != 0;
    }

    inline BaseSet BaseSet::set_difference(uint8_t a, uint8_t b) {
        return BaseSet(a & ~b);
    }

    inline BaseSet BaseSet::complement() {
        return BaseSet(0x0f_8 & ~set);
    }


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

        void emplace_read_column(std::vector<Read> &&col);

        const std::vector<std::vector<Read>> &get_read_columns() const;

        uint8_t get_ref() const;
    };

    typedef struct _Annotation{
        std::vector<int> genotype;
        std::vector<int> cnt_read;
        std::vector<int> cnt_tumor;

        _Annotation(int num_samples) : genotype(num_samples, 0), cnt_read(num_samples, 0), cnt_tumor(num_samples, 0) {}
    } Annotation;
    
}

#endif //MOSS_TYPES_H
