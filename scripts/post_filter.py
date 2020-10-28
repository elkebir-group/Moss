import numpy as np
from pysam import VariantFile
import argparse, sys

def apply_filters(ifile: str, ofile: str, threshold=100, normal_name="normal"):
    with VariantFile(ifile) as vcf_in:
        # Header
        header = vcf_in.header
        header.add_line('##FILTER=<ID=CLUSTER,Description="The site is clutered with at least 2 neighbors within 100bp">')
        with VariantFile(ofile, 'w', header=header) as vcf_out:
            for contig in vcf_in.header.contigs:
                buffer = []
                for rec in vcf_in.fetch(contig):
                    buffer.append(rec)
                    if len(buffer) == 3:
                        if 0 <= buffer[2].pos - buffer[0].pos <= threshold:
                            tumor_read_counts = np.array(
                                [[s["TCOUNT"] for k, s in buf.samples.items() if k!=normal_name] for buf in buffer])
                            is_clustered = np.any(np.all(tumor_read_counts>0, axis=0))
                            # print(contig, [buffer[i].pos for i in range(3)])
                            if is_clustered:
                                buffer[0].filter.add("CLUSTER")
                                buffer[1].filter.add("CLUSTER")
                                buffer[2].filter.add("CLUSTER")
                        vcf_out.write(buffer[0])
                        buffer.pop(0)
                else:
                    # end of for loop, write remaining reads
                    if len(buffer) > 0:
                        for rec in buffer:
                            vcf_out.write(rec)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input result VCF file")
    parser.add_argument("-o", "--output", help="output VCF file")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    apply_filters(args.input, args.output, normal_name="N1")
