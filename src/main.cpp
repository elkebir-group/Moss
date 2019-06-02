#include <iostream>
#include <getopt.h>
#include <vector>
#include <sys/stat.h>
#include "io/bam_io.h"
#include "core/SnvCaller.h"
#include "core/types.h"
#include <array>


const option long_options[] =
        {
                {"bam",     required_argument, nullptr,       'b'},
                {"ref",     required_argument, nullptr,       'r'},
                {"loci",    required_argument, nullptr,       'l'},
                {"tau",     required_argument, nullptr,       't'},
                {nullptr,   no_argument,       nullptr,       0}
        };

void print_help() {
    std::cout <<
        "\nMultiple Sample Somatic variant caller\n\n"
        "-h, --help             show this help message and exit\n"
        "-b, --bam <BAM>        a pair of original and realigned BAM files for\n"
        "                       one sample, can be addressed multiple times to\n"
        "                       specify\n"
        "-r, --ref <FASTA>      reference FASTA file\n"
        "-l, --loci <LOCI>      candidate loci files\n"
        "-t, --tau <TAU>        optional threshold for somatic score, default is -3"
        "--verbose              set verbose flag\n"
        "--brief                set brief flag, mutual exclusive to --verbose\n";
    exit(1);
}

bool is_file_exist(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}


int main(int argc, char **argv) {
    int c;
    const char* const short_opts = "b:r:l:t:h";


    // variables
    std::string ref_file;
    std::vector<std::string> bam_files;
    std::vector<std::string> bam_idx;
    std::vector<std::string> loci_files;
    double tau(-3.0);

    /* getopt_long stores the option index here. */
    int option_index = 0;
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

            case 'l':
                loci_files.emplace_back(std::string(optarg));
                break;

            case 't':
                tau = std::stod(optarg);
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

    if (!is_file_exist(ref_file)) {
        std::cerr << "Error: Failed to open reference file " << ref_file << std::endl;
        return 1;
    } else {
        std::cout << "Reference FASTA file: " << ref_file << std::endl;
    }
    std::cout << "BAM files:\t";
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
    std::cout << std::endl << "Loci files:\t";
    for (const auto &lociFile : loci_files) {
        std::cout << lociFile << '\t';
    }
    std::cout << "tau = " << tau << std::endl;

    // TODO: use htslib to parse files
    // TODO: merge loci files
    moss::SnvCaller caller(4);
    moss::BamStreamer streamer(ref_file, bam_files);
    std::array<int, 10> pos = {583, 631, 726, 760, 961, 990, 1270, 1507, 1705, 1743};
    auto p = pos[4];
    moss::Pileups col = streamer.get_column("demo20", p);
    auto array = col.get_read_columns();
    moss::BaseSet normal;
    uint8_t tumor;
    auto log_proba_non_soma = caller.calling(col, normal, tumor);
    std::cout<< "Pos: " << p << '\t' << "Prob: " << -10 * log_proba_non_soma << '\t' << seq_nt16_str[tumor] << std::endl;
    exit(0);
}