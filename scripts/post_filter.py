import numpy as np
from pysam import VariantFile
import argparse
import sys


def fallback(rec, normal_name):
    moss_support = sum(
        [s["TCOUNT"] > 0 for k, s in rec.samples.items() if k != normal_name])
    if moss_support > 1:
        if rec.info["NUMPASS"] > 0:
            if "LOW_QUAL" in rec.filter:
                rec.filter.clear()
                rec.filter.add("LOW_QUAL")
            else:
                rec.filter.clear()
                rec.filter.add("PASS")
    else:
        if rec.info["NUMPASS"] > 0:
            rec.filter.clear()
            rec.filter.add("PASS")
        else:
            rec.filter.clear()
            rec.filter.add("FALLBACK")


def apply_TIN(rec):
    if rec.info["TIN"] < 20:
        rec.filter.add("LOW_TIN")


def apply_filters(ifile: str, ofile: str, threshold=100, normal_name="normal", just_moss=False):
    with VariantFile(ifile) as vcf_in:
        # Header
        header = vcf_in.header
        header.add_line(
            '##FILTER=<ID=CLUSTER,Description="The site is clutered with at least 2 neighbors within 100bp">')
        header.add_line(
            '##FILTER=<ID=LOW_TIN,Description="TIN score is lower than 20">')
        header.add_line('##FILTER=<ID=FALLBACK,Description="Fallback to single-sample caller and NUMPASS=0">')
        header.add_line(f"##normal_sample={normal_name}")
        for sample in header.samples:
            if sample != normal_name:
                header.add_line(f"##tumor_sample={sample}")
        with VariantFile(ofile, 'w', header=header) as vcf_out:
            for contig in vcf_in.header.contigs:
                buffer = []
                for rec in vcf_in.fetch(contig):
                    buffer.append(rec)
                    if len(buffer) == 3:
                        if 0 <= buffer[2].pos - buffer[0].pos <= threshold:
                            tumor_read_counts = np.array(
                                [[s["TCOUNT"] for k, s in buf.samples.items() if k != normal_name] for buf in buffer])
                            is_clustered = np.any(np.all(tumor_read_counts > 0, axis=0))
                            # print(contig, [buffer[i].pos for i in range(3)])
                            if is_clustered:
                                buffer[0].filter.add("CLUSTER")
                                buffer[1].filter.add("CLUSTER")
                                buffer[2].filter.add("CLUSTER")
                        apply_TIN(buffer[0])
                        if not just_moss:
                            fallback(buffer[0], normal_name)
                        vcf_out.write(buffer[0])
                        buffer.pop(0)
                else:
                    # end of for loop, write remaining reads
                    if len(buffer) > 0:
                        for rec in buffer:
                            apply_TIN(rec)
                            if not just_moss:
                                fallback(rec, normal_name)
                            vcf_out.write(rec)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input result VCF file")
    parser.add_argument("-o", "--output", help="output VCF file")
    parser.add_argument("--normal-name", help="name of the normal sample")
    parser.add_argument("--just-moss", help="do not depend on single-sample caller result",
                        default=False, action='store_true')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    apply_filters(args.input, args.output,
                  normal_name=args.normal_name, just_moss=args.just_moss)
