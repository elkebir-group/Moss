//
// Created by Chuanyi Zhang on 2019-06-03.
//

#include <fstream>
#include <sstream>
#include <vector>
#include "loci.h"
#include "../io/vcf_io.h"

moss::MapContigLoci moss::merge_loci(std::vector<std::string> filenames) {
    MapContigLoci loci;
    for (const auto &filename : filenames) {
        std::string target,
                position;
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line, '\t')) {
            auto temp = std::stringstream(line);
            std::getline(temp, target, ':');
            std::getline(temp, position, ':');
            auto item = loci.find(target);
            if (item != loci.end()) {
                item->second.emplace(std::make_pair(
                    std::stoul(position),
                    Aggregate(false, 0)));
            } else {
                loci.insert(std::make_pair(target, std::map<locus_t, Aggregate>{{std::stoul(position), {false, 0}}}));
            }
        }
        file.close();
    }
    return loci;
}

moss::MapContigLoci moss::merge_vcf(std::vector<std::string> filenames) {
    MapContigLoci loci;
    for (const auto &filename : filenames) {
        VcfReader vcf(filename);
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

