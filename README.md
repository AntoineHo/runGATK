# runGATK
GATK variant calling pipeline using python

- [Environment and installation](#install)
- [Modules](#modules)
  - [Prepare module](#prepare)
  - [Trim module](#trim)
  - [Align module](#align)
  - [Call module](#call)
  - [Genotype module](#genotype)
- [Exemple usage: rungatk pipeline](#pipeline)
- [Exemple usage: manual run](#exemple)

## <a name="install"></a> 1. Environment and installation

### 1.1 Install conda and bioconda
See here: https://bioconda.github.io/user/install.html
- `conda` is required to install `numpy`, and `pyyaml`
- `bioconda` is required to install `bwa`, `sambamba`, `gatk4`, `picard`, `fastp`

### 1.2 Create environment & clone repository
```bash
conda/mamba create -n rungatk numpy pyyaml bwa sambamba gatk4 picard fastp
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

### <a name="trim"></a> 2.2 `trim`
This module trims PE reads libraries and removes adapters with `fastp`. In order it:
1. Reads an input tab separated file with 2 columns (respectively path to R1.fq.gz and path to R2.fq.gz). The reads can be (un)compressed.
2. Create `fastp` jobs and runs them in parallel

```bash
usage: runGATK.py trim [-h] [-t THREADS] [-p PROCESSES] [-fo FASTP_OPTIONS] Reads
positional arguments:
  Reads                <STRING> A file containing a list of filepaths to the reads to trim.
optional arguments:
  -h, --help           show this help message and exit
  -t, --threads        <INT> Maximum threads for fastp process. Default: [4]
  -p, --processes      <INT> Maximum parallel fastp processes in pool. Each process uses the number of threads defined by -t. Default: [4]
  -fo, --fastp-options <STRING> fastp filtering options. Default: ['--detect_adapter_for_pe --trim_poly_g --compression 9']
```

### <a name="align"></a> 2.3 `align`
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

### <a name="call"></a> 2.4 `call`
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

### <a name="genotype"></a> 2.5 `genotype`
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

## <a name="pipeline"></a> 3. Exemple usage: rungatk pipeline

The `pipeline` module creates a bash script file using a yaml formatted configuration file (see template below or download given file). This allows the user to simply input desired options and libraries, use the module and then run the bash script. This is done in one config file edition and 3 commands.

### 4.1 Preparing a configuration file

Exemple configuration file:
```yaml
# Dependencies must be in $PATH, use conda and bioconda for convenience
rungatk: /path/to/runGATK.py
reference: /path/to/reference.fasta
# Output directory (relative path from the directory where rungatk pipeline is called)
output_directory: output
# Options for trimming input reads
trim_options: # fastp trimming options
  skip_trimming: false # set to true to skip trimming
  njobs: 4 # Number of parallel trimming jobs (if num libraries < njobs then only num libraries in parallel)
  threads_per_job: 4 # Number of threads used by one job (total threads used = njobs * threads_per_job)
  fastp_options: '--detect_adapter_for_pe --trim_poly_g --compression 9'
align_options: # bwa & sambamba options
  threads: 20 # Alignement threads to use (one job at a time)
  ram_per_thread: 4 # RAM per thread for sambamba 
  mapq: 20 # minimum MAPQ to keep alignement (sambamba view)
prepare_options: # Options for chunking the input reference
  window_size: 50 # in Kbps
call_options: # GATK HaplotypeCaller options (check GATK4 manual for more information)
  njobs: 4 # Number of parallel HaplotypeCaller jobs
  ht: 4
  pim: 'NONE'
  java_options: '-Xmx4G' # Java Virtual Machine options (you may want to increase RAM here)
  he: 0.01 # heterozygosity
  ihe: 0.001 # indel heterozygosity
  mrpas: 50
  fi: ''
  erc: 'BP_RESOLUTION'
  om: 'EMIT_ALL_ACTIVE_SITES'
genotype_options: # GATK GenomicsDBImport and GenotypeGVCFs options (check GATK4 manual for more information)
  njobs: 4 # Number of parallel genotyping jobs (if num contigs < njobs then only num contigs in parallel) 
  sp: 4 # Sample in parallel for GenomicsDBImport
  java_options: '-Xmx4G'
  fi: ''
  he: 0.01
  ihe: 0.001
  not_allsites: false # Set this to 'true' if you are only looking for variant sites and not reference homozygous sites
samples:
  D4A3: # SAMPLE NAME that will be in the final VCF column header (do not use spaces)
    lib1: # PE library for sample D4A3. Note: the library name is not important for the pipeline (do not use spaces)
      R1: /path/to/input/D4A3.library.HiSeq2500.lane1.R1.fastq.gz # FWD
      R2: /path/to/input/D4A3.library.HiSeq2500.lane1.R2.fastq.gz # REVERSE
    lib2:
      R1: /path/to/input/D4A3.library.HiSeq2500.lane2.R1.fastq.gz
      R2: /path/to/input/D4A3.library.HiSeq2500.lane2.R2.fastq.gz 
  D5B3:
    lib1: # Only one library used for D5B3
      R1: /path/to/input/D5B3.lane1.R1.fastq.gz
      R2: /path/to/input/D5B3.lane1.R2.fastq.gz
```

### 4.2 Running the pipeline

```bash
conda activate rungatk
python runGATK.py pipeline config.yaml # The output directory is a relative path from the folder where this command launched
```

Move to the output directory where 2 files can be found :
- `pipeline.sh`
- `sample_reads.list`

It is important that you run pipeline.sh in this directory
```bash
cd /path/to/output
chmod +x pipeline.sh
./pipeline.sh
```

The final VCF file(s) are: `/path/to/output/jointgenotyping/merged.*.vcf`

## <a name="exemple"></a> 4. Exemple usage: manual run

### 4.1 Preparing reference

```bash
conda activate rungatk

wd=/path/to/workdir # path to working directory where runGATK was cloned
rungatk=$wd/runGATK/runGATK.py # path to runGATK.py

ref=/path/to/reference.fasta # path to reference file in .fasta format

mkdir -p $wd/input # create an input folder

cp $ref $wd/input # or ln -s $ref $wd/input # copy or symlink reference to the input folder

python $rungatk prepare -ws 500 $ref # Run prepare module
```

### 4.2 Aligning 2 libraries for 3 samples on the reference
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

### 4.3 Calling 3 samples
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

### 4.4 Joint genotyping of 3 samples
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
