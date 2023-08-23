# Assembly Free Linkage Analysis Pipeline

[![DOI](https://zenodo.org/badge/264294570.svg)](https://zenodo.org/badge/latestdoi/264294570)

## About

The Assembly Free Linkage Analysis Pipeline (AFLAP) is designed to build genetic maps using *k*-mers. This pipeline takes in raw reads, compares the composition (Jellyfish hashes) of the parents and identifies uniquely segregating *k*-mers. It has been tested on *Arabidopsis thaliana* and *Bremia lactucae* producing linkage groups coherent with independently generated genome assemblies. AFLAP may be applied to any organism for which a mapping population is available. We note that using low coverage sequence is less than ideal for an assembly-free approach as a large amount of noise is introduced due to missing data.

## Install

To install AFLAP, please clone this GitHub repo.

```bash
git clone https://github.com/kfletcher88/AFLAP.git
```

## Prerequisites

The following packages are required to run AFLAP:

- python v.3.10
- pandas v.1.5.3
- matplotlib v.3.7.1
- Jellyfish v.2.3.0
- ABySS v.2.3.7
- LepWrap v.4.0.1

All of the required dependencies can be easily installed by creating an environment via the provided yaml file (`AFLAP.yml`). It is highly recommended for the user to create the environment through conda.

All LepMap3 software is included through the LepWrap package.

## Custom Pedigree File

The pedigree file outlines the pedigree of the cross to be analyzed and the location of the reads which are going to be analyzed. This is done by mandatory tab delimited fields.

- Field 1 is the label which will be retained in all downstream analysis.
- Field 2 is the generation. This field will inform AFLAP what type of analysis is to be run:
  - parent (F0) = 0
  - child (F1) = 1
  - grandchild (F2) = 2.
- Field 3 is the location of the read file.

One individual may be split over multiple lines in instances where multiple read files are present. AFLAP will combine lines that share identical labels (Field 1).

### Additional Fields

For parents (F0), field 4 and 5 can be used to manually set the k-mer cutoffs used for marker assembly. If these fields are not set, AFLAP will try to estimate these cutoffs, which may not be perfect.\
Note AFLAP will plot the curve with the cut-offs supplied or calculated. These can be used to edit the Pedigree file and rerun AFLAP. Rerunning AFLAP is efficient as it will reuse all applicable, previously calculated results.\
If the parents are low coverage or GBS then the estimates may be very far off.\
For progeny (F1 and F2), field 4 and 5 can be used to specify the parents. Not currently implemented, but is foundational for a potential multi-cross analysis enhancement.

An example pedigree is available in:

```bash
./example/Pedigree.txt
```

It looks like this:

```bash
#label  F       Read    LowerCoverage/Parent1   UpperCoverage/Parent2
Col     0       ReadsEBI/SRR5882797_1.fastq.gz  32      105
Col     0       ReadsEBI/SRR5882797_2.fastq.gz  32      105
Ler     0       ReadsEBI/SRR3166543_1.fastq.gz  40      177
Ler     0       ReadsEBI/SRR3166543_2.fastq.gz  40      177
ERR1432507      2       ReadsEBI/ERR1432507_1.fastq.gz  Col     Ler
ERR1432507      2       ReadsEBI/ERR1432507_2.fastq.gz  Col     Ler
ERR1432492      2       ReadsEBI/ERR1432492_1.fastq.gz  Col     Ler
ERR1432492      2       ReadsEBI/ERR1432492_2.fastq.gz  Col     Ler
```

Currently AFLAP is launched using the Python script `AFLAP.py` which accepts multiple options.

```text
$ python3 AFLAP.py <-P> ...

Options:
  -P Pedigree file (required). See AFLAP README for more
     information.
  -m K-mer size. Default [31].
  -t Threads for JELLYFISH counting. Default [4].
  -r Individual to remove. All other options will be ignored.
  -L LOD score - Will run LepMap3 with minimum LOD. Default
     [2].
  -d Lower boundary for marker cut off. Can be used to filter
     for segregation distortion. Default [0.2].
  -D Upper boundary for marker cut off. Can be used to filter
     for segregation distortion. Default [0.8].
  -f Limit for how many XX can exist in a row. If surpassed
     then sequence is not considered for analysis. Default
     [None].
  -x Run with low coverage parameters.
  -n Minimum number of linkage groups needed to continue F2
     LepMap3 analysis. Default [10].
  -U Maximum number of markers to output in the genotype
     tables output under ./AFLAP_Results/
```

## Intermediate Results

All temporary files are stored in AFLAP_tmp. This can be pretty sizeable depending on the biology of the organism understudy, especially if you are running different parameters. Once you have generated a genetic map, then I recommend you delete this directory. When AFLAP is initiated it will look for this and it will not regenerate files already present, so storing this data can save time.

## Final Results

There are multiple points which AFLAP can be stopped:

1. **I just want a genotype table to export to my favorite mapping software.** This is generated in step 5, so users can terminate the software once this stage is completes if desired. The genotype table will be saved in AFLAP_Results, with the suffix GT.tsv.
2. **I want a LepMap3 compatible genotype table that I will use to run LepMap3 independently.** This is generated in step 6 and is the default stopping position of AFLAP. It will be saved in AFLAP_Results, with the suffix ForLepMap3.tsv.
3. **I want LepMap3 results.** By specifying minimal LOD scores when running `AFLAP.py`, using `-l` AFLAP will run LepMap3, provided it is found in the ThirdParty directory. In running this AFLAP will calculate the marker order for every linkage group which contains more than 1% of the total number of markers. In addition, AFLAP will process the LepMap3 results to produce a final file detailing the Marker ID, Linkage group assigned, cM position and marker sequence. This can be used to generate a marker file to map to a genome assembly very quickly. Marker sequences can be used to compare between runs. MarkerIDs can be used for comparison, provided the same AFLAP_tmp directory is used (e.g. in case the user wants to change the progeny).\
The LepMap3 results are located in AFLAP_Results/LOD#, depending on the minimum LOD score provided for F1 analysis progeny or the LOD score that identifies the use-specified number of linkage groups (argument `-n`) for F2 analysis.\
The Map and marker file is located in AFLAP_Results suffixed LOD#.txt.

## Frequently Asked Questions

Q: My run crashed halfway through, do I have to start again?\
A: Yes, but AFLAP will detect intermediate files and not overwrite them, therefore it should progress to the last point quite quickly.

Q: I want to run multiple times with different coverage cutoffs/k-mer sizes, is it possible?\
A: Yes, AFLAP stores parameters in the file names and will write new files if new parameters are detected.\
For example JELLYFISH results calculated at 21 or 31 base pairs will be suffixed jf21 or jf31 respectively. The resulting genotype tables will have "m21" or "m31" in there file names.

Q: I wish to add individuals to my Pedigree file, do I have to start in a new directory/run the full pipeline on every individual?\
A: No, AFLAP will be able to use old results for previously generated data and generate the required files for the new sequences. Just add these to the Pedigree file.

Q: I want to exclude individuals, should I delete intermediate files?\
A: No, just provide a Pedigree file without those individuals. The genotype table is directed with the Pedigree file, so will only build a table for progeny indicated with in.

Q: I have added sequences to an individual and thus added lines to the pedigree file, will AFLAP detect this?\
A: No, though I might add something in the future. In the meantime, you can supply `-r` to `AFLAP.py` to remove intermediate files for specific progeny individuals. Rerunning AFLAP will then automatically recalculate those intermediate files. E.G:

```bash
AFLAP.py -P Pedigree.txt -m 31 -t 8
AFLAP.py -P Pedigree.txt -r ProgenyInd1
```
