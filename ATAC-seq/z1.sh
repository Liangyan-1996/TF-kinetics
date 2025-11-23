thread="8"

for sample in $(cat zzsample.txt);do
sbatch --job-name=job --partition=cu --nodes=1 --ntasks-per-node=$thread --output=%j.out --error=%j.err --wrap="
sh /public/home/liunangroup/liangyan/pipeline/mypipeline/ATAC-seq/ATAC-seq.sh $sample $thread mm
"
done
