fastp="/public/home/liunangroup/liangyan/software/miniconda3/envs/QC/bin/fastp"
STAR="/public/home/liunangroup/liangyan/software/miniconda3/envs/rnaseq/bin/STAR"
sambamba="/public/home//liunangroup/liangyan/software/miniconda3/envs/chipseq/bin/sambamba"
STAR_index="/public/home/liunangroup/liangyan/Genome/ensembl/hg38/STAR_full"
GTF="/public/home/liunangroup/liangyan/Genome/ensembl/hg38/Homo_sapiens.GRCh38.110.chr.gtf"

thread="32"

mkdir QC
mkdir mapping

for sample in `cat zzsample.txt`;do
sbatch --job-name=job --partition=cu --nodes=1 --ntasks-per-node=$thread --output=%j.out --error=%j.err --wrap="
$fastp \
   -i raw/${sample}_1.fastq.gz \
   -I raw/${sample}_2.fastq.gz \
   -o QC/${sample}_1.fq.gz \
   -O QC/${sample}_2.fq.gz \
   -w $thread

$STAR \
    --genomeDir $STAR_index \
    --runThreadN $thread \
    --outBAMsortingThreadN $thread \
    --outFilterMismatchNmax 20 \
    --outFilterMatchNminOverLread 0.4 \
    --outFilterScoreMinOverLread 0.4 \
    --readFilesIn QC/${sample}_1.fq.gz QC/${sample}_2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes nM MD NH \
    --outFileNamePrefix mapping/${sample}. \
    --quantMode GeneCounts

$sambamba index -t $thread mapping/${sample}.Aligned.sortedByCoord.out.bam

rm -r QC/${sample}_*.gz
"
done
