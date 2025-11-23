#!/bin/bash
sample="$1"
thread="$2"

fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
bowtie2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bowtie2"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/samtools"
sambamba="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/sambamba"
bedtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bedtools"
macs2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/macs2"
bamCoverage="/public/home/liunangroup/liangyan/software/miniconda3/bin/bamCoverage"
computeMatrix="/public/home/liunangroup/liangyan/software/miniconda3/bin/computeMatrix"
plotHeatmap="/public/home/liunangroup/liangyan/software/miniconda3/bin/plotHeatmap"
python="/public/home/liunangroup/liangyan/software/miniconda3/bin/python"
ScriptDir="/public/home/liunangroup/liangyan/pipeline/mypipeline/script"
homer="/public/home/liunangroup/liangyan/software/homer/bin/findMotifsGenome.pl"

bowtie2index="/public/home/liunangroup/liangyan/Genome/gencode/hg38/bowtie2_index/genome"
blacklist="/public/home/liunangroup/liangyan/Genome/blacklist/hg38.blacklist.bed"
GenomeSize="2913022398"

workdir=`pwd`
cd $workdir
mkdir qc mapping log peak_narrow peak_broad motif

$fastp -i raw/${sample}_1.fastq.gz -I raw/${sample}_2.fastq.gz -o qc/${sample}_1.fq.gz -O qc/${sample}_2.fq.gz -w $thread
$bowtie2 --end-to-end --very-sensitive -x $bowtie2index -p $thread -X 700 -I 10 --no-mixed --no-discordant --no-unal -1 qc/${sample}_1.fq.gz -2 qc/${sample}_2.fq.gz 2> $workdir/log/$sample.log |\
$sambamba view -S -F "ref_name != 'chrM'" -f bam -t $thread -o mapping/$sample.align.bam /dev/stdin

cd mapping
# sorting
$sambamba sort -t $thread -m 10G -o $sample.sorted.bam $sample.align.bam
# remove PCR duplicates
$sambamba markdup -r -t $thread $sample.sorted.bam $sample.dedup.bam 2>> $workdir/log/$sample.log
# remove blacklist reads
$samtools view -@ $thread -o $sample.blacklist.bam -U $sample.rmblacklist.bam -L $blacklist $sample.dedup.bam
# final sorting
$sambamba sort -t $thread -m 10G -o $sample.bam $sample.rmblacklist.bam

# coverage
$bamCoverage -b $sample.bam -o $sample.bw -p $thread --normalizeUsing RPGC --effectiveGenomeSize $GenomeSize

# peak calling
$macs2 callpeak -t $sample.bam -f BAMPE --gsize hs -n $sample -q 0.01 --keep-dup all --outdir $workdir/peak_narrow
$macs2 callpeak -t $sample.bam -f BAMPE --gsize hs -n $sample -q 0.01 --keep-dup all --broad --broad-cutoff 0.05 --outdir $workdir/peak_broad

# heatmap
$computeMatrix reference-point -S $sample.bw -R $workdir/peak_narrow/${sample}_summits.bed -o $sample.gz -a 3000 -b 3000 -p $thread
$plotHeatmap --matrixFile $sample.gz -out $sample.png --dpi 300 --plotFileFormat png --yMin 0

# motif enrichment
$homer $workdir/peak_narrow/${sample}_summits.bed hg38 $workdir/motif/$sample -size -100,100 -p $thread

# postqc
# length distribution
$sambamba sort -t $thread -m 10G -N -o $sample.sorted_by_name.bam $sample.bam
$bedtools bamtobed -bedpe -i $sample.sorted_by_name.bam > $sample.bedpe
$python $ScriptDir/CUTTag/frag_len_distribution.py $sample.bedpe $sample.fragment_length.png

# MT percent
line_align=$($sambamba view -t $thread -c $sample.sorted.bam)
line_MT=$($sambamba view -t $thread -c $sample.sorted.bam chrM)
MT_pct=$(echo -e "scale=5;$line_MT / $line_align * 100"|bc)
echo -e "MT proportion: $MT_pct%" >> $workdir/log/$sample.log

# FRiP
line_clean=$($sambamba view -t $thread -c $sample.bam)
line_inPeak=$($sambamba view -t $thread -c -L $workdir/peak_narrow/${sample}_peaks.narrowPeak $sample.bam)
FRiP=$(echo -e "scale=5;$line_inPeak / $line_clean"|bc)
echo -e "Clean Fragments: $line_clean" >> $workdir/log/$sample.log
echo -e "FRiP: $FRiP" >> $workdir/log/$sample.log

# remove tmp
rm $sample.align.bam* $sample.sorted.bam* $sample.dedup.bam* $sample.blacklist.bam* $sample.rmblacklist.bam* $sample.sorted_by_name.bam*

# merge results
cd $workdir
sh /public/home/liunangroup/liangyan/pipeline/mypipeline/script/CUTTag/generate_report.sh
