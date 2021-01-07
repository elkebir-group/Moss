import yaml
import argparse, sys, os
import time


def main(args):
    with open(args.config, 'r') as config_file, open(args.output, 'w') as bash_file:
        folder = os.path.dirname(args.config) + '/'
        config = yaml.safe_load(config_file)
        ref = folder + config["ref"]
        vcfs = ['-i ' + folder + v for v in config["vcfs"]]
        caller = config["caller"]
        candidate = folder + config["candidate"]
        normal_name = config["normal_name"]
        output = folder + config["output"]
        output_gz = output + ".gz"
        output_filtered = output[:output.rfind('.')] + ".post_filter.vcf"
        moss_log = folder + f"moss.{int(time.time())}.log"
        if config["realigned_bams"]:
            if len(config["bams"]) != len(config["realigned_bams"]):
                raise ValueError(f"error in config file `{args.config}`, length of `bams` ({len(config['bams'])}) is different than `realigned_bams` ({len(config['realigned_bams'])})")
            list_of_bams = []
            for bam, realign in zip(config["bams"], config["realigned_bams"]):
                if realign:
                    list_of_bams.append('-b ' + folder + bam + ' -R ' + folder + realign)
                else:
                    list_of_bams.append('-b ' + folder + bam)
        else:
            list_of_bams = ['-b ' + folder + bam for bam in config["bams"]]

        commands = []
        commands.append("#!/bin/bash\n")
        commands.append(f"python /moss_scripts/union_candidates.py -t {caller} {' '.join(vcfs)} --normal-name {normal_name} -o {candidate}\n")
        commands.append(f"moss -r {ref} {' '.join(list_of_bams)} -l {candidate} -m 4 -t -0.693 --ignore0 --grid-size 200 -o {output} > {moss_log}\n")
        commands.append(f"bgzip {output}\n")
        commands.append(f"tabix {output_gz}\n")
        commands.append(f"python /moss_scripts/post_filter.py --normal-name {normal_name} -i {output_gz} -o {output_filtered}\n")
        bash_file.writelines(commands)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Run Moss")
    parser.add_argument("-c", "--config", help="config file containing list of bam files and other info")
    parser.add_argument("-o", "--output", help="output path for bash script")
    args = parser.parse_args(None if sys.argv[1:] else ["-h"])

    main(args)
