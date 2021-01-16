# runGATK
GATK variant calling pipeline using python

- [Environment and installation](#install)
- [Modules](#modules)
  - [Prepare module](#prepare)
  - [Align module](#align)
  - [Call module](#call)
  - [Genotype module](#genotype)
- [Exemple usage](#exemple)

## <a name="install"></a> 1. Environment and installation

### 1.1 Install conda and bioconda
See here: https://bioconda.github.io/user/install.html
- `conda` is required to install `numpy` and `biopython`
- `bioconda` is required to install `bwa`, `sambamba`, `samtools`, `gatk4`, `picard`
Note: it is possible that in the future some of these requirements will drop (e.g.: `numpy`, `biopython`) and others will be added (e.g.:`fastp`)

### 1.2 Create environment & clone repository
```bash
conda create -n rungatk biopython numpy bwa sambamba samtools picard gatk4
conda activate rungatk
cd $workdir # Folder of your project
git clone https://github.com/AntoineHo/runGATK.git
```

## <a name="modules"></a> 2. Module usage

### <a name="prepare"></a> 2.1 `prepare`
This module creates window intervals for parallel calling with `HaplotypeCaller` and genotyping with `GenotypeGVCF`

```bash
usage: runGATK.py prepare [-h] [-ws WINDOW_SIZE] Reference
positional arguments:
  Reference           <STRING> A fasta file containing the reference genome.
optional arguments:
  -h, --help          show this help message and exit
  -ws, --window-size  <INT> Window size (in kb) to split genome. Default: [100]
```

### <a name="align"></a> 2.2 `align`
This module creates a .bam file from the input reads. In order it:
1. Aligns the reads library per library and lane by lane (with automatic read-groups naming) with `bwa mem`
2. Filters bad quality alignment and converts SAM to BAM file with `sambamba view`
3. Sorts temporary .bam files with `sambamba sort`
4. Merges all temporary filtered .bam files with `sambamba merge`
5. Mark duplicates in merged filtered .bam file with `sambamba markdup`
6. Sort final merged, filtered and marked .bam file with `sambamba sort`
7. Remove all temporary .sam and .bam files

```bash
usage: runGATK.py align [-h] [-t THREADS] [-r RAM] [-m MAPQ] -kp Sample Reference R1 R2
positional arguments:
  Sample            <STRING> sample name for output and adding to read group. Default: 'sample1'
  Reference         <STRING> A path to the reference genome fasta file.
  R1                <STRING> A comma separated list of FastQ files FORWARD (R1) (can be compressed).
  R2                <STRING> A comma separated list of FastQ files REVERSE (R2) (can be compressed).
optional arguments:
  -h, --help        show this help message and exit
  -t, --threads     <INT> Maximum threads to use. Default: [4]
  -r, --ram         <INT> Number of GB of RAM per thread to use. Default: [4]
  -m, --mapq        <INT> Sambamba view filtering to keep only alignements with MAPQ >= <INT>. Default: [20]
  -kp, --keep-temp  Do not remove intermediate steps alignments (cannot be before positional argument). Default: False
```

### <a name="call"></a> 2.3 `call`
This module creates .g.vcf files from the sample bam files. In order it:
1. Calls samples .bam files in parallel using intervals defined by the `prepare` module with `GATK HaplotypeCaller`
2. Merge intervals .g.vcf into a single merged.g.vcf file per sample with `picard MergeVcfs`
3. Remove temporary files.

```bash
usage: runGATK.py call [-h] [-pim {NONE,CONSERVATIVE,HOSTILE,AGGRESSIVE}] [-jo JAVA_OPTIONS] [-he HETEROZYGOSITY] [-p PROCESSES] [-ht HMM_THREADS] [-mrpas MAX_READS_PER_ALIGNMENT_START] [-ihe INDEL_HETEROZYGOSITY] [-fi FOUNDER_ID] [-erc {BP_RESOLUTION,GVCF}] [-om {EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_ACTIVE_SITES}] [-kp] [-dr] [-bo] Sample Reference

positional arguments:
  Sample                                  <STRING> sample name for output and SAME AS "runGATK.py align"
  Reference                               <STRING> A path to the reference genome fasta file
optional arguments:
  -h, --help                              show this help message and exit
  -pim, --pcr-indel-model                 <STRING> Argument to pass to HaplotypeCaller --pcr-indel-model. Valid options: (NONE,CONSERVATIVE,HOSTILE,AGGRESSIVE). Default: ['NONE']
  -jo, --java-options                     <STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: ['-Xmx4G']
  -he, --heterozygosity                   <FLOAT> Heterozygosity value to pass to HaplotypeCaller. Default: [0.01]
  -p, --processes                         <INT> Maximum processes in pool to use. Default: [4]
  -ht, --hmm-threads                      <INT> Number of threads to use for HMM. Default: [4] (per process).
  -mrpas, --max-reads-per-alignment-start <INT> See GATK doc. Default: [50] (set to 0 to deactivate).
  -ihe, --indel-heterozygosity            <FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: [0.0001]
  -fi, --founder-id                       <STRING> Founder population sample ID. Default: empty
  -erc, --emit-ref-confidence             <STRING> See HaplotypeCaller documentation Emit Ref Confidence. Valid options: (BP_RESOLUTION, GVCF). Default: ['BP_RESOLUTION']
  -om, --output-mode                      <STRING> See HaplotypeCaller documentation output mode. Valid options: (EMIT_VARIANTS_ONLY,EMIT_ALL_CONFIDENT_SITES,EMIT_ALL_ACTIVE_SITES). Default: ['EMIT_ALL_ACTIVE_SITES']
  -kp, --keep-temp                        Do not remove intermediate steps files (cannot be before positional argument). Default: False
  -dr, --dry-run                          Only parse arguments for testing and debugging commands. Default: False
```

### <a name="genotype"></a> 2.4 `genotype`
This module creates a final .vcf file from the samples merged.g.vcf files. In order it:
1. Imports .g.vcf files from all the samples to a database with `GATK GenomicsDBImport` <- best option
1. (Alternatively) Combines .g.vcf into a single multi-sample .g.vcf with `GATK CombineGVCFs` <- much longer runtime
2. Genotype multisample g.vcf files in parallel chromosome (contig) per chromosome with `GATK GenotypeGVCFs`. Warning: by default it creates 2 files (raw.vcf and allsites.vcf). If you are only interested in the sites that have variants, use the `--not-all` flag (faster).
3. Merge genotyped chromosome .vcf files into a single final .vcf file containing all samples and all sites (or variant only).
4. Remove temporary .vcf files.

```bash
usage: runGATK.py genotype [-h] [-p PROCESSES] [-sp SAMPLE_PROCESSES] [-bs BATCH_SIZE] [-jo JAVA_OPTIONS] [-fi FOUNDER_ID] [-he HETEROZYGOSITY] [-ihe INDEL_HETEROZYGOSITY] [-na] [-kp] [--use-combine] Reference Samples Output

positional arguments:
  Reference                     <STRING> A path to the output directory of "runGATK.py prepare"
  Samples                       <STRING> A comma separated list of samples (directory must contain call/merged.g.vcf)
  Output                        <STRING> An output directory path for the joint calling.
optional arguments:
  -h, --help                    show this help message and exit
  -p, --processes               <INT> Maximum threads to use. Default: [4]
  -sp, --sample-processes       <INT> Maximum threads to read samples simultaneously (GenomicsDBImport only). Default: [4]
  -bs, --batch-size             <INT> Default batch size for GenomicsDBImport. Default: [50]
  -jo, --java-options           <STRING> Java Virtual Machine options (Ram Per Process is defined here). Default: ['-Xmx4G']
  -fi, --founder-id             <STRING> Founder population sample ID. Default: "" (empty)
  -he, --heterozygosity         <FLOAT> Heterozygosity value to pass to HaplotypeCaller. Default: [0.01]
  -ihe, --indel-heterozygosity  <FLOAT> Indel heterozygosity value to pass to HaplotypeCaller. Default: [0.0001]
  -na, --not-all                Do not use allsites option in GenotypeGVCFs (in addition to normal vcf output). Default: will output allsites vcf and normal vcf output.
  -kp, --keep-temp              Do not remove intermediate steps files (cannot be before positional argument). Default: False
  --use-combine                 Use CombineGVCFS instead of GenomicsDBImport
```

## <a name="exemple"></a> 3. Exemple usage

### 3.1 Preparing reference

```bash
conda activate rungatk

wd=/path/to/workdir # path to working directory where runGATK was cloned
rungatk=$wd/runGATK/runGATK.py # path to runGATK.py

ref=/path/to/reference.fasta # path to reference file in .fasta format

mkdir -p $wd/input # create an input folder

cp $ref $wd/input # or ln -s $ref $wd/input # copy or symlink reference to the input folder

python $rungatk prepare -ws 500 $ref # Run prepare module
```

### 3.2 Aligning 2 libraries for 3 samples on the reference
Note: reads pre-processing is not implemented (yet?) in the pipeline. I usually run `fastQC` and `fastp` before using the reads.

```bash
conda activate rungatk

wd=/path/to/workdir # path to working directory where runGATK was cloned
rungatk=$wd/runGATK/runGATK.py # path to runGATK.py

ref=/path/to/reference.fasta # path to reference file in .fasta format

cd $wd # output will be stored in $wd/$sample for each run
for sample in sample1 sample2 sample3; do
L1R1=/path/to/$sample/library1.R1.fq.gz
L1R2=/path/to/$sample/library1.R2.fq.gz
L2R1=/path/to/$sample/library2.R1.fq.gz
L2R2=/path/to/$sample/library2.R2.fq.gz

# Give in the same order the R1 and then the R2 libs (comma separated):
python $rungatk align $sample $ref $L1R1,$L2R1 $L1R2,$L2R2
done
```

### 3.3 Calling 3 samples
```bash
conda activate rungatk

wd=/path/to/workdir # path to working directory where runGATK was cloned
rungatk=$wd/runGATK/runGATK.py # path to runGATK.py

ref=/path/to/reference.fasta # path to reference file in .fasta format
processes=12 # Run 12 jobs in parallel

cd $wd # Output will be stored in $wd/$sample/call for each sample
for sample in sample1 sample2 sample3; do
python $rungatk call -p $processes $sample $ref # Run HaplotypeCaller in on multiple intervals in parallel
done
```

### 3.4 Joint genotyping of 3 samples
```bash
conda activate rungatk

wd=/path/to/workdir # path to working directory where runGATK was cloned
rungatk=$wd/runGATK/runGATK.py # path to runGATK.py

ref=/path/to/reference.fasta # path to reference file in .fasta format

processes=12 # Run 12 jobs max in parallel (if less than 12 contigs in reference then max jobs in parallel = nbr of contigs else max 12)
cd $wd
# Note: directories $wd/sample1 $wd/sample2 and $wd/sample3 must exist and contain a subdirectory "call" with a merged .g.vcf file in it
# Note: output directory is chosen and is not especially the working directory for this command
python $rungatk genotype -p $processes $ref sample1,sample2,sample3 $wd/output
```
