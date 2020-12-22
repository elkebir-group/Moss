#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 chuanyi5 <chuanyi5@illinois.edu>
#
# Distributed under terms of the MIT license.

"""
Generate candidates VCF from Strelka2, Mutect2 tumor VCFs
"""

import allel
import argparse
import numpy as np
import sys


def strelka2(inputs) -> dict:
    # chrom: pos: normal_genotype
    chrom_pos_gt = dict()
    cnt = 0
    for ifile in inputs:
        in_vcf = allel.read_vcf(ifile, fields='*')
        for chrom, pos, ref, alt, nt, is_pass in zip(in_vcf["variants/CHROM"], in_vcf["variants/POS"], in_vcf["variants/REF"], in_vcf["variants/ALT"], in_vcf["variants/NT"], in_vcf["variants/FILTER_PASS"]):
            chrom = str(chrom)
            pos = int(pos)
            alt = alt[0]
            num_pass = int(is_pass)

            if nt == 'ref':
                normal = ref + ref
            elif nt == 'het':
                normal = ref + alt
            elif nt == 'hom':
                normal = alt + alt
            else:
                continue

            if chrom in chrom_pos_gt:
                if pos in chrom_pos_gt[chrom]:
                    chrom_pos_gt[chrom][pos]["num_pass"] += num_pass
                    if chrom_pos_gt[chrom][pos]["gt"][0] == chrom_pos_gt[chrom][pos]["gt"][1] and normal[0] != normal[1]:
                        chrom_pos_gt[chrom][pos]["gt"] = normal
                        cnt += 1
                else:
                    chrom_pos_gt[chrom][pos] = {"gt": normal, "num_pass": num_pass}
            else:
                chrom_pos_gt[chrom] = {pos: {"gt": normal, "num_pass": num_pass}}
    print(f"Disagreement on normal: {cnt} times.")
    return chrom_pos_gt


def mutect2(inputs, normal_name) -> dict:
    # chrom: pos: normal_genotype
    chrom_pos_gt = dict()
    cnt = 0
    cnt_het_hom = 0
    for ifile in inputs:
        # ["variants/CHROM", "variants/POS", "variants/REF", "variants/ALT", "calldata/GT"]
        in_vcf = allel.read_vcf(ifile, fields='*')
        idx_normal = np.argwhere(in_vcf["samples"] == normal_name)[0][0]
        zipped = zip(in_vcf["variants/CHROM"][in_vcf["variants/is_snp"]],
            in_vcf["variants/POS"][in_vcf["variants/is_snp"]],
            in_vcf["variants/REF"][in_vcf["variants/is_snp"]],
            in_vcf["variants/ALT"][in_vcf["variants/is_snp"]],
            in_vcf["calldata/GT"][in_vcf["variants/is_snp"]],
            in_vcf["variants/FILTER_PASS"][in_vcf["variants/is_snp"]],
            in_vcf["variants/FILTER_artifact_in_normal"][in_vcf["variants/is_snp"]])
        for chrom, pos, ref, alt, gt, is_pass, is_artifact in zipped:
            if is_artifact:
                continue
            chrom = str(chrom)
            pos = int(pos)
            num_pass = int(is_pass)
            alt = alt[0]
            ref_alt = ref + alt
            normal = ref_alt[gt[idx_normal][0]] + ref_alt[gt[idx_normal][1]]
            if gt[idx_normal][0] != 0 or gt[idx_normal][1] != 0:
                cnt_het_hom += 1

            if chrom in chrom_pos_gt:
                if pos in chrom_pos_gt[chrom]:
                    chrom_pos_gt[chrom][pos]["num_pass"] += num_pass
                    if chrom_pos_gt[chrom][pos]["gt"][0] == chrom_pos_gt[chrom][pos]["gt"][1] and normal[0] != normal[1]:
                        cnt += 1
                        chrom_pos_gt[chrom][pos]["gt"] = normal
                else:
                    chrom_pos_gt[chrom][pos] = {"gt": normal, "num_pass": num_pass}
            else:
                chrom_pos_gt[chrom] = {pos: {"gt": normal, "num_pass": num_pass}}
    print(f"Disagreement on normal: {cnt} times.")
    print(f"Not ref: {cnt_het_hom} times.")
    return chrom_pos_gt


def write_vcf(dict_chrom_pos_gt: dict, output, input_files, is_split: bool):
    header = """##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=NUMPASS,Number=1,Type=Integer,Description=\"Number of samples that pass the base caller\">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\n"""
    if is_split:
        for chrom, pos_info in dict_chrom_pos_gt.items():
            prefix = 'chr' if not chrom.startswith("chr") else ''
            out_filename = output[:output.rfind('.')] + f".{prefix}{chrom}.vcf"
            with open(out_filename, "w") as ofile:
                ofile.write(header)
                chunk = []
                cnt = 0
                for pos, info in pos_info.items():
                    gt = info["gt"]
                    num_pass = info["num_pass"]
                    info = f"NUMPASS={num_pass}"
                    if gt[0] == gt[1]:
                        chunk.append(
                            f"{chrom}\t{pos}\t.\t{gt[0]}\t.\t.\t.\t{info}\tGT\t0|0\n")
                    cnt += 1
                    if cnt % 1000 == 0:
                        ofile.writelines(chunk)
                        chunk = []
                        cnt = 0
                ofile.writelines(chunk)
    else:
        with open(output, "w") as ofile:
            ofile.write(header)
            for chrom, pos_info in dict_chrom_pos_gt.items():
                chunk = []
                cnt = 0
                for pos, info in pos_info.items():
                    gt = info["gt"]
                    num_pass = info["num_pass"]
                    info = f"NUMPASS={num_pass}"
                    if gt[0] == gt[1]:
                        chunk.append(
                            f"{chrom}\t{pos}\t.\t{gt[0]}\t.\t.\t.\t{info}\tGT\t0|0\n")
                    cnt += 1
                    if cnt % 1000 == 0:
                        ofile.writelines(chunk)
                        chunk = []
                        cnt = 0
                ofile.writelines(chunk)

def main(args):
    if args.input_files is not None:
        with open(args.input_files) as ifile:
            if args.input is None:
                args.input = []
            for line in ifile:
                args.input.append(line.strip())
    
    
    if args.tool[0].lower() == 'm':
        if args.normal_name is None:
            exit()
        chrom_pos_gt = mutect2(args.input, args.normal_name)
    elif args.tool[0].lower() == 's':
        chrom_pos_gt = strelka2(args.input)
    write_vcf(chrom_pos_gt, args.output, args.input_files, args.split)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input tumor result VCF file", action="append")
    parser.add_argument("--normal-name", help="name of the normal sample in the VCF file, only used for Mutect")
    parser.add_argument("-f", "--input-files", help="input tumor result VCF file list")
    parser.add_argument("-t", "--tool", help="[M|m|Mutect] or [S|s|Strelka]", required=True)
    parser.add_argument("-o", "--output", help="output candidates VCF file")
    parser.add_argument("--split", help="split output VCF by chromosomes", action="store_true")
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
