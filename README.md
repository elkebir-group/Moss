# Moss

Moss is a multi-sample somatic single nucleotide variant (SNV) calling tool aiming for discovering variants with a low variant allele frequency that repeatedly appears in several samples.
Moss works as an extension to existing single-sample somatic variant callers and improves the recall meanwhile maintains high precision.
Moss takes as input the BAM files of multiple samples and corresponding VCF output of the single-sample caller.

![Figure](doc/overview.png)

## Contents

  1. [Compilation instructions](#compilation)
     * [Dependencies](#dep)
     * [Compilation](#comp)
  2. [Usage instructions](#usage)

<a name="compilation"></a>

## Compilation instructions

<a name="dep"></a>

### Dependencies

Moss is written in C++11 and thus requires a modern C++ compiler (GCC >= 4.8.1, or Clang). In addition, Moss has the following dependencies.

* [HTSlib](https://github.com/samtools/htslib/releases) (>=1.7)
* [CMake](http://www.cmake.org) (>= 3.9)
* Python (>= 3.6)
* [Scikit-allel](https://pypi.org/project/scikit-allel/)

<a name="comp"></a>

### Compilation

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
```

If `HTSlib` is not in the system path, CMake may not be able to find it. Users then need to manually set the path for `htslib` using `ccmake`:

```bash
ccmake ..
```

Then finally make Moss.

```bash
make
```

The compilation results in the executable `moss`.

<a name="usage"></a>

## Usage instructions

Moss works on top of other somatic variant calling methods, such as Strelka2 and Mutect2.
We assume you have already run the base variant caller as their manual suggested and get the VCF files of each sample.
Then run the python script `scripts/union_candidates.py` to generate a VCF file of candidates loci as an input to Moss. For example:

```bash
python scripts/union_candidates.py -f <list_of_VCF.list> --normal-name <NORMAL> -t Mutect -o <output.vcf>
```

To run Moss, you need a reference genome FASTA file, BAM files for normal and tumor samples, realigned BAM files (optional but recommended), and a candidate loci VCF.
For example, after you've built `moss` in the `build/` directory, you can run the toy example in `data/`:

``` bash
./moss -r ../data/demo20.fa -b ../data/normal.sort.bam -R ../data/empty.bam -b ../data/clone0.spike.sort.bam -R ../data/empty.bam -b ../data/clone1.spike.sort.bam -R ../data/empty.bam -b ../data/clone2.spike.sort.bam -R ../data/empty.bam -b ../data/clone3.spike.sort.bam -R ../data/empty.bam -l ../data/candidates.chrdemo20.vcf -m 4 -t -0.693 --ignore0 --grid-size 200 -o example.vcf
```
