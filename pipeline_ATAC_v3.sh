#!/bin/bash
set -euo pipefail

#NOTE: samtools >=1.16
###VARIABLES
#run
#FQ1=/mnt/data/rawdata/g99/R2100710240008/FT150041020/L01/FT150041020_L01_read_1.fq.gz
#FQ2=/mnt/data/rawdata/g99/R2100710240008/FT150041020/L01/FT150041020_L01_read_2.fq.gz
#OUTDIR=/mnt/nvme/home/chris/OneStep_Francesco_072026_v2
THREADS=10
#ref
#GENOMEDIR=/mnt/nvme/home/chris/resources/genomes/Human/hg38.primary_assembly.fa
GENOMEDIR=/mnt/nvme/home/chris/resources/genomes/Human/hg38.primary_assembly
#BARCODELIST=/mnt/data/user/chris/ATAC/OneStep_Francesco_072026_v2/samplesheet.csv
BARCODELIST=/mnt/nvme/home/chris/resources/barcodes/smartATAC_ATAC_bcs_v1_mgikit_percell.csv
#tools
SAMTOOLS=/mnt/nvme/home/chris/miniconda3/envs/samtools/bin/samtools
#BWA=/mnt/nvme/home/chris/programs/bwa-mem2-2.2.1_x64-linux/bwa-mem2
#BWA=/mnt/nvme/home/chris/programs/bwa/bwa
BOWTIE=/mnt/nvme/home/chris/programs/bowtie2-2.5.5-linux-x86_64/bowtie2
MGIKIT=/mnt/nvme/home/chris/programs/mgikit-2.1.1/bins/mgikit
FASTP=/mnt/nvme/home/chris/programs/fastp
#SEQKIT=/mnt/nvme/home/chris/programs/seqkit
SINTO=/mnt/nvme/home/chris/miniconda3/envs/sinto/bin/sinto
#LOCAL=no

# --- Required (no default) ---
FQ1=""
FQ2=""
OUTDIR=""

usage() {
    cat <<EOF
Usage: $0 --fq1 FILE --fq2 FILE --outdir DIR [options]

Required:
  --fq1 FILE            Path to read 1 FASTQ
  --fq2 FILE            Path to read 2 FASTQ
  --outdir DIR          Output directory

Optional (default shown):
  --threads N           Number of threads (default: $THREADS)
  --genomedir FILE      Reference genome fasta (default: $GENOMEDIR)
  --barcodelist FILE    Barcode list csv (default: $BARCODELIST)
  --samtools PATH       Path to samtools binary (default: $SAMTOOLS)
  --bowtie PATH         Path to bowtie2 binary (default: $BOWTIE)
  --mgikit PATH         Path to mgikit binary (default: $MGIKIT)
  --fastp PATH          Path to fastp binary (default: $FASTP)
  --sinto PATH          Path to sinto binary (default: $SINTO)
  -h, --help            Show this help message
EOF
    exit 1
}
#--seqkit PATH         Path to seqkit binary (default: $SEQKIT)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fq1)         FQ1="$2"; shift 2 ;;
        --fq2)         FQ2="$2"; shift 2 ;;
        --outdir)      OUTDIR="$2"; shift 2 ;;
        --threads)     THREADS="$2"; shift 2 ;;
        --genomedir)   GENOMEDIR="$2"; shift 2 ;;
        --barcodelist) BARCODELIST="$2"; shift 2 ;;
        --samtools)    SAMTOOLS="$2"; shift 2 ;;
        --bowtie)      BOWTIE="$2"; shift 2 ;;
        --mgikit)      MGIKIT="$2"; shift 2 ;;
        --fastp)       FASTP="$2"; shift 2 ;;
        #--seqkit)      SEQKIT="$2"; shift 2 ;;
        --sinto)       SINTO="$2"; shift 2 ;;
       #--local)       LOCAL="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

# --- Validate required args ---
if [[ -z "$FQ1" || -z "$FQ2" || -z "$OUTDIR" ]]; then
    echo "Error: --fq1, --fq2, and --outdir are required." >&2
    usage
fi

# --- Echo config ---
echo "Running scATAC-seq pipeline with following parameters:"
echo "FQ1:         $FQ1"
echo "FQ2:         $FQ2"
echo "OUTDIR:      $OUTDIR"
echo "THREADS:     $THREADS"
echo "GENOMEDIR:   $GENOMEDIR"
echo "BARCODELIST: $BARCODELIST"
echo "SAMTOOLS:    $SAMTOOLS"
echo "BOWTIE:      $BOWTIE"
echo "MGIKIT:      $MGIKIT"
echo "FASTP:       $FASTP"
#echo "SEQKIT:      $SEQKIT"
echo "SINTO:       $SINTO"
#echo "LOCAL:        $LOCAL"

START=$(date +%s)
get_timing () {
    END=$(date +%s)
    ELAPSED=$((END - START))
    printf "Step done. Time elapsed: %dm %ds\n" $((ELAPSED / 60)) $((ELAPSED % 60))
}

if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

###STEP1: Barcode handling/MGIKIT
echo "STEP1: Barcode handling/MGIKIT"
$MGIKIT demultiplex --read1 $FQ1 --read2 $FQ2 --sample-sheet $BARCODELIST --output $OUTDIR/demux --template i720 --run 01 --threads $THREADS --instrument R2100400190055
get_timing

###STEP2: renaming fastq headers/SEQKIT
echo "STEP2: renaming fastq headers"
nextcall () {
    i=$1
    if [ -e $OUTDIR/demux/"$i"_*_R1_001.fastq.gz ]; then
        #$SEQKIT replace -p '^(\S+) (.+):([^:]+)$' -r '${3}:${1} ${2}' --threads 1 -o $OUTDIR/demux/"$i"_R1_001.renamed.fastq.gz $OUTDIR/demux/"$i"_*_R1_001.fastq.gz
        #$SEQKIT replace -p '^(\S+) (.+):([^:]+)$' -r '${3}:${1} ${2}' --threads 1 -o $OUTDIR/demux/"$i"_R2_001.renamed.fastq.gz $OUTDIR/demux/"$i"_*_R2_001.fastq.gz
        zcat $OUTDIR/demux/"$i"_*_R1_001.fastq.gz | awk -v var="$i" 'NR%4==1{ sub(/^@/, "@" var ":") } {print}' | gzip > $OUTDIR/demux/"$i"_R1_001.renamed.fastq.gz
        zcat $OUTDIR/demux/"$i"_*_R2_001.fastq.gz | awk -v var="$i" 'NR%4==1{ sub(/^@/, "@" var ":") } {print}' | gzip > $OUTDIR/demux/"$i"_R2_001.renamed.fastq.gz
    fi
}

#cut -f1 $BARCODELIST | grep -v 'sample_id' > $OUTDIR/sample_list.txt
find $OUTDIR/demux/ -size +100k -type f -name '*R1_001.fastq.gz' -printf "%f\n" | sed 's/_/\t/g' | awk '{print $1"_"$2}' | grep -v Undetermined > $OUTDIR/sample_list.txt
for i in $(cat $OUTDIR/sample_list.txt); do
    #njobs=$(pgrep -cu "$(whoami)" -f "seqkit" || true)
    njobs=$(pgrep -cu "$(whoami)" -f "awk" || true)
    while [ "$njobs" -ge "$THREADS" ]; do
        sleep 10
        njobs=$(pgrep -cu "$(whoami)" -f "seqkit" || true)
    done
    nextcall "$i" &
done
wait
cat $OUTDIR/demux/*R1_001.renamed.fastq.gz > $OUTDIR/demux/ATAC.renamed.read1.fastq.gz
cat $OUTDIR/demux/*R2_001.renamed.fastq.gz > $OUTDIR/demux/ATAC.renamed.read2.fastq.gz
get_timing

###STEP3: Trimming/FASTP
echo "STEP3: Trimming/FASTP"
if [ ! -d "$OUTDIR/trimmed" ]; then
    mkdir $OUTDIR/trimmed
fi
$FASTP \
    --in1 $OUTDIR/demux/ATAC.renamed.read1.fastq.gz \
    --in2 $OUTDIR/demux/ATAC.renamed.read2.fastq.gz \
    --out1 $OUTDIR/trimmed/ATAC.trimmed.read1.fastq.gz \
    --out2 $OUTDIR/trimmed/ATAC.trimmed.read2.fastq.gz \
    --trim_front1 27 --trim_front2 19 --trim_poly_x \
    --cut_tail --cut_right_window_size 10 --qualified_quality_phred 20 \
    --json $OUTDIR/trimmed/ATAC.fastp.json --html $OUTDIR/trimmed/ATAC.fastp.html --thread $THREADS

rm $OUTDIR/demux/*.fastq.gz
get_timing


###STEP4: Alignmend/BWA + mate fixing/SAMOOLS + cell barcode tagging/SINTO + read sorting/SAMTOOLS
echo "STEP4: Alignmend/BOWTIE + read sorting/SAMTOOLS"
if [ ! -d "$OUTDIR/aligned" ]; then
    mkdir $OUTDIR/aligned
fi
#if [ "$LOCAL" == "yes" ]; then
    $BOWTIE --local --threads $THREADS -X 2000 -x $GENOMEDIR -1 $OUTDIR/trimmed/ATAC.trimmed.read1.fastq.gz -2 $OUTDIR/trimmed/ATAC.trimmed.read2.fastq.gz | $SAMTOOLS fixmate -m -u - -  | $SINTO nametotag --tag CB --bam - | $SAMTOOLS sort -@ $THREADS -m 2G -o $OUTDIR/aligned/ATAC.aligned.sorted.tmp.bam
#else
#    $BOWTIE --threads $THREADS -X 2000 -x $GENOMEDIR -1 $OUTDIR/trimmed/ATAC.trimmed.read1.fastq.gz -2 $OUTDIR/trimmed/ATAC.trimmed.read2.fastq.gz | $SAMTOOLS fixmate -m -u - -  | $SINTO nametotag --tag CB --bam - | $SAMTOOLS sort -@ $THREADS -m 2G -o $OUTDIR/aligned/ATAC.aligned.sorted.tmp.bam
#fi
get_timing


###STEP5: Cell-type-aware duplicate marking/SAMTOOLS
echo "STEP5: Cell-type-aware duplicate marking/SAMTOOLS"
$SAMTOOLS markdup \
    -@ $THREADS \
    --barcode-tag CB \
    -S -s \
    $OUTDIR/aligned/ATAC.aligned.sorted.tmp.bam \
    $OUTDIR/aligned/ATAC.aligned.sorted.bam \
    2> $OUTDIR/aligned/ATAC.dup_metrics.txt

rm $OUTDIR/aligned/ATAC.aligned.sorted.tmp.bam
get_timing


###STEP6: index/SAMTOOLS
echo "STEP6: index bam file/SAMTOOLS"
$SAMTOOLS index \
    -@ "$THREADS" \
    $OUTDIR/aligned/ATAC.aligned.sorted.bam
get_timing

###STEP7: fragment file/SINTO
echo "STEP7: create fragment file/SINTO"
$SINTO fragments --use_chrom "^[1-9|M|chr]" --bam $OUTDIR/aligned/ATAC.aligned.sorted.bam --nproc $THREADS --barcodetag CB --fragments $OUTDIR/ATAC.frags.bed 
sort -k 1,1 -k2,2n $OUTDIR/ATAC.frags.bed > $OUTDIR/ATAC.frags.sorted.bed
/mnt/nvme/home/chris/miniconda3/envs/samtools/bin/bgzip $OUTDIR/ATAC.frags.sorted.bed
/mnt/nvme/home/chris/miniconda3/envs/samtools/bin/tabix -p bed $OUTDIR/ATAC.frags.sorted.bed.gz
get_timing

echo "DONE"
