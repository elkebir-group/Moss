//
// Created by Chuanyi Zhang on 2019-06-03.
//

#include <fstream>
#include <sstream>
#include <vector>
#include "loci.h"
#include "../io/vcf_io.h"

moss::MapContigLoci moss::import_loci(std::string filename) {
    MapContigLoci loci;
    VcfReader<LociRecData> vcf(filename);
    for (auto &&chrom_pos : vcf.get_records()) {
        auto item = loci.find(chrom_pos.first);
        if (item != loci.end()) {
            for (auto &&pos_rec : chrom_pos.second) {
                auto found_pos = item->second.find(pos_rec.first);
                if (found_pos != item->second.end()) {
                    found_pos->second.is_pass = pos_rec.second.num_pass > 0;
                    found_pos->second.num_pass = pos_rec.second.num_pass;
                } else {
                    item->second.emplace(std::make_pair(
                        pos_rec.first,
                        Aggregate(pos_rec.second.num_pass > 0, pos_rec.second.num_pass)));
                }
            }
        } else {
            loci.insert(std::make_pair(chrom_pos.first, std::map<locus_t, Aggregate>()));
            for (auto &&pos_rec : chrom_pos.second) {
                loci.at(chrom_pos.first).emplace(std::make_pair(
                    pos_rec.first,
                    Aggregate(pos_rec.second.num_pass > 0, pos_rec.second.num_pass)));
            }
        }
    }
    return loci;
}

moss::MapContigLoci moss::merge_vcf(std::vector<std::string> filenames) {
    MapContigLoci loci;
    for (const auto &filename : filenames) {
        VcfReader<RecData> vcf(filename);
        for (auto &&chrom_pos : vcf.get_records()) {
            auto item = loci.find(chrom_pos.first);
            if (item != loci.end()) {
                for (auto &&pos_rec : chrom_pos.second) {
                    auto found_pos = item->second.find(pos_rec.first);
                    if (found_pos != item->second.end()) {
                        found_pos->second.is_pass = pos_rec.second.is_pass | found_pos->second.is_pass;
                        found_pos->second.num_pass += static_cast<unsigned>(pos_rec.second.is_pass);
                    } else {
                        item->second.emplace(std::make_pair(
                            pos_rec.first,
                            Aggregate(pos_rec.second.is_pass, static_cast<unsigned>(pos_rec.second.is_pass))));
                    }
                }
            } else {
                loci.insert(std::make_pair(chrom_pos.first, std::map<locus_t, Aggregate>()));
                for (auto &&pos_rec : chrom_pos.second) {
                    loci.at(chrom_pos.first).emplace(std::make_pair(
                        pos_rec.first,
                        Aggregate(pos_rec.second.is_pass, static_cast<unsigned>(pos_rec.second.is_pass))));
                }
            }
        }
    }
    return loci;
}

