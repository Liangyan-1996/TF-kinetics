sample="$1"
thread="$2"
specie="$3"

fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
bowtie2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bowtie2"
samtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/samtools"
sambamba="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/sambamba"
bedtools="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/bedtools"
macs2="/public/home/liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/macs2"
alignmentSieve="/public/home/liunangroup/liangyan/software/miniconda3/bin/alignmentSieve"
computeMatrix="/public/home/liunangroup/liangyan/software/miniconda3/bin/computeMatrix"
plotHeatmap="/public/home/liunangroup/liangyan/software/miniconda3/bin/plotHeatmap"
bamCoverage="/public/home/liunangroup/liangyan/software/miniconda3/bin/bamCoverage"
Rscript="/public/home/liunangroup/liangyan/software/miniconda3/envs/ATAC/bin/Rscript"
python="/public/home/liunangroup/liangyan/software/miniconda3/bin/python"
ScriptDir="/public/home/liunangroup/liangyan/pipeline/mypipeline/script"

if [ "$specie" = "hs" ]; then
    bowtie2index="/public/home/liunangroup/liangyan/Genome/gencode/hg38/bowtie2_index/genome"
    blacklist="/public/home/liunangroup/liangyan/Genome/blacklist/hg38.blacklist.bed"
    TSS="/public/home/liunangroup/liangyan/Genome/gencode/hg38/hg38.TSS.bed"
    GenomeSize="2913022398"
elif [ "$specie" = "mm" ]; then
    bowtie2index="/public/home/liunangroup/liangyan/Genome/gencode/mm10/bowtie2_index/genome"
    blacklist="/public/home/liunangroup/liangyan/Genome/blacklist/mm10.blacklist.bed"
    TSS="/public/home/liunangroup/liangyan/Genome/gencode/mm10/mm10.TSS.bed"
    GenomeSize="2652783500"
else
    echo "error: unsupported species $specie"
    exit 1
fi

workdir=$(pwd)
cd $workdir
mkdir qc mapping log peak

if [[ -e raw/${sample}_1.fastq.gz ]] && [[ -e raw/${sample}_2.fastq.gz ]];then
library="2"
echo "PE data"

$fastp -i raw/${sample}_1.fastq.gz -I raw/${sample}_2.fastq.gz -o qc/${sample}_1.fq.gz -O qc/${sample}_2.fq.gz -w $thread
$bowtie2 --end-to-end --very-sensitive -x $bowtie2index -p $thread -X 2000 --no-mixed --no-discordant --no-unal -1 qc/${sample}_1.fq.gz -2 qc/${sample}_2.fq.gz 2>> $workdir/log/$sample.log |\
$sambamba view -S -f bam -t $thread -o mapping/$sample.align.bam /dev/stdin

elif [[ -e raw/$sample.fastq.gz ]];then
library="1"
echo "SE data"

$fastp -i raw/$sample.fastq.gz -o qc/$sample.fq.gz -w $thread -A
$bowtie2 --end-to-end --very-sensitive -x $bowtie2index -p $thread --no-unal -U qc/$sample.fq.gz 2>> $workdir/log/$sample.log |\
$sambamba view -S -f bam -t $thread -o mapping/$sample.align.bam /dev/stdin

else
echo 'error, no raw data'

fi

cd mapping
# sorting
$sambamba sort -t $thread -m 10G -o $sample.sorted.bam $sample.align.bam
# remove MT
$sambamba view -F "ref_name != 'chrM'" -f bam -t $thread -o $sample.rmMT.bam $sample.sorted.bam
# remove PCR duplicate
$sambamba markdup -r -t $thread $sample.rmMT.bam $sample.dedup.bam 2>> $workdir/log/$sample.log
# remove blacklist
$samtools view -@ $thread -o $sample.blacklist.bam -U $sample.rmblacklist.bam -L $blacklist $sample.dedup.bam
# final sorting
$sambamba sort -t $thread -m 10G -o $sample.bam $sample.rmblacklist.bam
# coverage
$bamCoverage -b $sample.bam -o $sample.bw -p $thread --normalizeUsing RPGC --effectiveGenomeSize $GenomeSize --extendReads 100

# peak calling
if [ "$library" == "2" ];then
$macs2 callpeak -t $sample.bam -f BAMPE --gsize $specie -n $sample -q 0.01 --keep-dup all --outdir $workdir/peak
elif [ "$library" == "1" ];then
$macs2 callpeak -t $sample.bam -f BAM --gsize $specie -n $sample -q 0.01  --shift -70 --extsize 140 --nomodel --keep-dup all --outdir $workdir/peak
fi

# TSS heatmap
$computeMatrix reference-point -S $sample.bw -R $TSS -a 3000 -b 3000 -p $thread -o $sample.gz
$plotHeatmap -m $sample.gz -out $sample.png --dpi 300 --plotFileFormat png --yMin 0

# MT percent
line_align=$($sambamba view -t $thread -c $sample.sorted.bam)
line_MT=$($sambamba view -t $thread -c $sample.sorted.bam chrM)
MT_pct=$(echo -e "scale=5;$line_MT / $line_align * 100"|bc)
echo -e "MT proportion: $MT_pct%" >> $workdir/log/$sample.log

# TSS enrichment score
$Rscript $ScriptDir/ATAC-seq/TSS_enrichment_score.R $sample.bam $specie >> $workdir/log/$sample.log

# length distribution
if [ "$library" == "2" ];then
$sambamba sort -t $thread -m 10G -N -o $sample.sorted_by_name.bam $sample.bam
$bedtools bamtobed -bedpe -i $sample.sorted_by_name.bam > $sample.bedpe
$python $ScriptDir/ATAC-seq/frag_len_distribution.py $sample.bedpe $sample.fragment_length.png
fi

# FRiP
line_clean=$($sambamba view -t $thread -c $sample.bam)
line_inPeak=$($sambamba view -t $thread -c -L $workdir/peak/${sample}_peaks.narrowPeak $sample.bam)
FRiP=$(echo -e "scale=5;$line_inPeak / $line_clean" | bc)
echo -e "Clean Fragments: $line_clean" >> $workdir/log/$sample.log
echo -e "FRiP: $FRiP" >> $workdir/log/$sample.log

# remove tmp
rm $sample.align.bam* $sample.sorted.bam* $sample.rmMT.bam* $sample.dedup.bam* $sample.blacklist.bam* $sample.rmblacklist.bam* $sample.sorted_by_name.bam*

# summary
cd $workdir
sh $ScriptDir/ATAC-seq/ATAC-seq_summarize.sh
