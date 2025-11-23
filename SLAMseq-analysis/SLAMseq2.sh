mkdir grandslam

thread="96"
mem=$(echo "$thread * 10" | bc)G
ls mapping/*.bam | sort > grandslam/zzsample.bamlist

sbatch --job-name=job --partition=fat --nodes=1 --ntasks-per-node=$thread --mem-per-cpu=10G --output=%j.out --error=%j.err --wrap="
gedi \
    -mem $mem \
    -e Slam \
    -genomic ensembl_hg38 \
    -nthreads $thread \
    -mode Unique \
    -highmem \
    -trim5p 10 \
    -trim3p 10 \
    -prefix grandslam/1.GrandSlam \
    -reads grandslam/zzsample.bamlist \
    -full \
    -plot \
    -D
"
