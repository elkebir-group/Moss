#include <iostream>
#include <getopt.h>
#include <vector>
#include <sys/stat.h>
#include <bitset>
#include "io/bam_io.h"
#include "io/loci.h"
#include "core/calling.h"
#include "core/types.h"
#include <array>
#include <cmath>
#include <chrono>
#include <iomanip>

struct Flags {
    int dry = 0;
    int filter_total = 0;
    int filter_vaf = 0;
    int ignore0 = 0;
} flags;

const option long_options[] =
    {
        {"bam",                 required_argument, nullptr,   'b'},
        {"ref",                 required_argument, nullptr,   'r'},
        {"output",              required_argument, nullptr,   'o'},
        {"normal",              required_argument, nullptr,   'n'},
        {"loci",                required_argument, nullptr,   'l'},
        {"vcf",                 required_argument, nullptr,   'v'},
        {"tau",                 required_argument, nullptr,   't'},
        {"mu",                  required_argument, nullptr,   'm'},
        {"max_dep",             required_argument, nullptr,   'd'},
        {"min-base-qual",       required_argument, nullptr,   'B'},
        {"min-mapping-qual",    required_argument, nullptr,   'M'},
        {"grid-size",           required_argument, nullptr,   'g'},
        {"filter-total",         no_argument,       &flags.filter_total, 1},
        {"filter-vaf",           no_argument,       &flags.filter_vaf, 1},
        {"ignore0",             no_argument,       &flags.ignore0, 1},
        {"dry",                 no_argument,       &flags.dry,  1},
        {nullptr,               no_argument,       nullptr,   0}
    };

void print_help() {
    std::cout <<
              "Multiple Sample Somatic variant caller\n"
              "\n"
              "Usage:\n"
              "  moss [options] -r <reference> -b <normal BAM> -b <tumor BAM>\n"
              "       -n <germline VCF> -v <candidate VCF> -o <output VCF>\n"
              "\n"
              "Compulsary arguments:\n"
              "  -b, --bam <BAM>        a pair of original and realigned BAM files for\n"
              "                         one sample, can be addressed multiple times\n"
              "  -r, --ref <FASTA>      reference FASTA file\n"
              "  -o, --output <OUT>     output VCF file\n"
              "  -n, --normal <NORMAL>  normal sample's germline VCF result\n"
              "  -v, --vcf <VCF>        tumor samples' somatic VCF result\n"
              "  -l, --loci <LOCI>      candidate loci files, 0-based\n"
              "\n"
              "General options\n"
              "  -h, --help             show this help message and exit\n"
              "\n"
              "Other options:\n"
              "  -t, --tau <TAU>        optional threshold for somatic score, default is 0\n"
              "  -m, --mu <MU>          set somatic mutation prior, 1-5x10^(-m), default is 1-5e-4\n"
              "  -B, --min-base-qual <MINBASE>\n"
              "                         minimum threshold for base quality, default is 13\n"
              "  -M, --min-mapping-qual <MINMAP>\n"
              "                         minimum threshold for mapping quality, default is 30\n"
              "  -d, --max_dep <depth>  max depth of the dataset, default is 500\n"
              "  -g, --grid-size <grid> number of grid points in [0, 1] for approximate integration,\n"
              "                         default is 20, must be non-negative\n"
              "  --dry                  set to dry run\n"
              "  --filter-total          set to filter variants with total tumor depth < 150\n"
              "  --filter-vaf            set to filter variants with VAF in any samples < 0.1\n"
              "  --ignore0              ignore samples with ref allele in all reads\n";
    exit(1);
}

bool is_file_exist(const std::string &name) {
    struct stat buffer;
    return (stat(name.c_str(), &buffer) == 0);
}


int main(int argc, char **argv) {
    int c;
    const char *const short_opts = "b:r:o:n:l:v:t:m:B:M:d:h";


    // variables
    std::string ref_file;
    std::string normal_vcf;
    std::string out_vcf;
    std::vector<std::string> bam_files;
    std::vector<std::string> bam_idx;
    std::vector<std::string> loci_files;
    std::vector<std::string> tumor_vcfs;
    float tau{0};
    double mu{1 - 5e-4};
    int max_depth{500},
        min_base_qual{13},
        min_mapping_qual{30},
        grid_size{20};


    /* getopt_long stores the option index here. */
    int option_index = 0;
    if (argc == 1) {
        print_help();
    }
    while (true) {
        c = getopt_long(argc, argv, short_opts, long_options, &option_index);
        if (c == -1) {
            break;
        }
        switch (c) {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != nullptr)
                    break;
                printf("option %s", long_options[option_index].name);
                if (optarg)
                    printf(" with arg %s", optarg);
                printf("\n");
                break;

            case 'b':
                bam_files.emplace_back(std::string(optarg));
                break;

            case 'r':
                ref_file = std::string(optarg);
                break;

            case 'o':
                out_vcf = std::string(optarg);
                break;

            case 'n':
                normal_vcf = std::string(optarg);
                break;

            case 'l':
                loci_files.emplace_back(std::string(optarg));
                break;

            case 'v':
                tumor_vcfs.emplace_back(std::string(optarg));
                break;

            case 't':
                tau = std::stof(optarg);
                break;

            case 'm':
                mu = 1 - 5 * pow(10, -std::stod(optarg));
                break;

            case 'B':
                min_base_qual = std::stoi(optarg);
                break;

            case 'M':
                min_mapping_qual = std::stoi(optarg);
                break;

            case 'g':
                grid_size = std::stoi(optarg);
                break;

            case 'd':
                max_depth = std::stoi(optarg);
                break;

            case '?':
            case 'h':
            default:
                /* getopt_long already printed an error message. */
                print_help();
        }
    }

    /* Print any remaining command line arguments (not options). */
    if (optind < argc) {
        printf("non-option ARGV-elements: ");
        while (optind < argc)
            printf("%s ", argv[optind++]);
        putchar('\n');
    }

    if (flags.dry == 1) {
        std::cout << "# Dry run!" << std::endl;
    }

    if (!is_file_exist(ref_file)) {
        std::cerr << "Error: Failed to open reference file " << ref_file << std::endl;
        return 1;
    } else {
        std::cout << "# Reference FASTA file: " << ref_file << std::endl;
    }
    std::cout << "# BAM files:\t";
    for (const std::string &bamFile : bam_files) {
        std::string indexFile = bamFile + ".bai";
        if (!is_file_exist(bamFile)) {
            std::cerr << "Error: Failed to open BAM file " << bamFile << std::endl;
            return 1;
        } else if (!is_file_exist(indexFile)) {
            std::cerr << "Error: Failed to open index file " << indexFile << std::endl;
            return 1;
        }
        std::cout << bamFile << '\t';
    }
    std::cout << std::endl << "# Loci files:\t";
    for (const auto &lociFile : loci_files) {
        if (!is_file_exist(lociFile)) {
            std::cerr << "Error: Failed to open loci file " << lociFile << std::endl;
            return 1;
        } else {
            std::cout << lociFile << '\t';
        }
    }
    std::cout << std::endl << "# Input VCF files:\t";
    if (!normal_vcf.empty()) {
        if (!is_file_exist(normal_vcf)) {
            std::cerr << "Error: Failed to open VCF file " << normal_vcf << std::endl;
            return 1;
        } else {
            std::cout << normal_vcf << '\t';
        }
    }
    for (const auto &tumor_vcf : tumor_vcfs) {
        if (!tumor_vcf.empty()) {
            if (!is_file_exist(tumor_vcf)) {
                std::cerr << "Error: Failed to open VCF file " << tumor_vcf << std::endl;
                return 1;
            } else {
                std::cout << tumor_vcf << '\t';
            }
        }
    }
    std::cout << std::endl << "# Output VCF file:\t";
    if (out_vcf.empty()) {
        std::cerr << std::endl << "Error: No VCF file specified. Use `-o output.vcf` option to specify." << std::endl;
        return 1;
    } else {
        std::cout << out_vcf << std::endl;
    }
    std::cout << "tau = " << tau << " mu = " << mu << std::endl;

    unsigned long num_samples = bam_files.size();
    unsigned long num_tumor_samples = bam_files.size() - 1;
//    auto loci = moss::merge_loci(loci_files);
    moss::MapContigLoci loci;
    if (tumor_vcfs.size() > 0) {
        loci = moss::merge_vcf(tumor_vcfs);
    } else if (loci_files.size() > 0) {
        loci = moss::merge_loci(loci_files);
    } else {
        std::cerr << std::endl
                  << "Error: No candidate VCF or loci file specified. Use `-v candidate.vcf` option or `-l candidate.loci` to specify."
                  << std::endl;
        return 1;
    }
    std::cout << "# Loci merged" << std::endl;
    auto start = std::chrono::system_clock::now();
    moss::SnvCaller caller(num_tumor_samples, normal_vcf, flags.ignore0, mu, max_depth, grid_size+1);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::chrono::duration<double> calling_time{0};
    std::cout << "# SnvCaller elapsed time: " << elapsed_seconds.count() << std::endl;
    std::cout << "## Chrom\tPos \t Prob \t Alt \t Genotype:TumorCount:Coverage" << std::endl;
    moss::Annotation anno(num_samples);
    start = std::chrono::system_clock::now();
    moss::VcfWriter writer{out_vcf, loci, num_tumor_samples, ref_file, bam_files, flags.filter_total, flags.filter_vaf, tau};
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    std::cout << "# Writer elapsed time: " << elapsed_seconds.count() << std::endl;
    auto now = std::time(nullptr);
    std::cout << "# Start time: " << std::put_time(std::localtime(&now), "%FT%T%z") << std::endl;
    for (const auto &chrom : loci) {
        start = std::chrono::system_clock::now();
        moss::BamStreamer streamer(ref_file, bam_files, moss::MapContigLoci{{chrom.first, chrom.second}}, min_base_qual, min_mapping_qual);
        end = std::chrono::system_clock::now();
        elapsed_seconds = end - start;
        std::cout << "# Streamer elapsed time: " << elapsed_seconds.count() << std::endl;
        for (const auto &l : chrom.second) {
            moss::Pileups col = streamer.get_column();
            if (flags.dry == 0) {
                // TODO: baseset
                start = std::chrono::system_clock::now();
                auto log_proba_non_soma = caller.calling(chrom.first, l, col, anno);
                end = std::chrono::system_clock::now();
                calling_time += end - start;

                writer.write_record(chrom.first, l, col.get_ref(), anno, num_samples);

                if (log_proba_non_soma < tau) {
                    std::cout << chrom.first << '\t' << l + 1 << '\t' << -10 * log_proba_non_soma << '\t'
                              << seq_nt16_str[anno.tumor_gt] << '\t';
                    for (size_t i = 0; i < num_samples; i++) {
                        std::cout << "0|" << anno.genotype[i] << ':'
                                  << anno.cnt_tumor[i] << ','
                                  << anno.cnt_read[i] << '\t';
                    }
                    std::cout << std::endl;
                }

            }
        }
    }
    std::cout << "# Calling elapsed time: " << calling_time.count() << std::endl;
    now = std::time(nullptr);
    std::cout << "# End time: " << std::put_time(std::localtime(&now), "%FT%T%z") << std::endl;
    return 0;
}

