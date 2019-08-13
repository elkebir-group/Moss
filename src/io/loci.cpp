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
                item->second.insert(std::stoul(position));
            } else {
                loci.insert(std::make_pair(target, std::set<unsigned long>({std::stoul(position)})));
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
                    item->second.insert(pos_rec.first);
                }
            } else {
                loci.insert(std::make_pair(chrom_pos.first, std::set<unsigned long>()));
                for (auto &&pos_rec : chrom_pos.second) {
                    loci.at(chrom_pos.first).insert(pos_rec.first);
                }
            }
        }
    }
    return loci;
}

