# Moss

A multi-sample somatic SNV caller

Usage:
``` bash
moss [options] -r <reference> -b <normal BAM> -b <tumor BAM>
       -n <germline VCF> -v <candidate VCF> -o <output VCF>
```

Moss works on top of other somatic variant calling methods, such as Strelka2 and Mutect2.
We assume you have already run the base variant caller as their manual suggested and get the VCF files of each sample.
Then with the python converting script located at `scripts/use_normal_gt.py`, you can generate the germline VCF file from the VCFs produced by the base caller.

```bash
python scripts/use_normal_gt.py -f <list_of_VCF.list> -t Mutect -o <output.vcf>
```


