#!/bin/bash

# 基因组预处理模块
# 包含基因组消化和Bowtie2索引构建

prepare_genome() {
    local genome_file="$1"
    local restriction_site="$2"
    local outdir="$3"
    local cpu="$4"
    
    blue "开始基因组预处理..."
    
    # 创建基因组目录
    local genome_dir="$outdir/genome"
    mkdir -p "$genome_dir"
    
    # 步骤1: 基因组消化
    blue "步骤1: 基因组消化"
    local fragments_file="$genome_dir/$(basename "${genome_file%.*}")_fragments.bed"
    
    if [[ ! -f "$fragments_file" ]]; then
        blue "使用digest_genome.py进行基因组消化..."
        python3 d:\download\下载\hic-2.1.0\bin\digest_genome.py \
            "$genome_file" \
            -r "$restriction_site" \
            -o "$fragments_file" \
            > "$outdir/logs/digest_genome.log" 2>&1
        
        if [[ $? -eq 0 ]]; then
            green "基因组消化完成: $fragments_file"
        else
            red "基因组消化失败,请检查日志: $outdir/logs/digest_genome.log"
            exit 1
        fi
    else
        yellow "基因组消化文件已存在,跳过: $fragments_file"
    fi
    
    # 步骤2: Bowtie2索引构建
    blue "步骤2: Bowtie2索引构建"
    local genome_name=$(basename "${genome_file%.*}")
    local bt2_index_base="$genome_dir/bowtie2_index/$genome_name"
    
    # 检查索引是否存在
    local index_files=("$bt2_index_base.1.bt2" "$bt2_index_base.2.bt2" "$bt2_index_base.3.bt2" "$bt2_index_base.4.bt2" "$bt2_index_base.rev.1.bt2" "$bt2_index_base.rev.2.bt2")
    local index_exists=true
    
    for index_file in "${index_files[@]}"; do
        if [[ ! -f "$index_file" ]]; then
            index_exists=false
            break
        fi
    done
    
    if [[ "$index_exists" == false ]]; then
        blue "构建Bowtie2索引..."
        mkdir -p "$genome_dir/bowtie2_index"
        
        bowtie2-build \
            --threads "$cpu" \
            "$genome_file" \
            "$bt2_index_base" \
            > "$outdir/logs/bowtie2_build.log" 2>&1
        
        if [[ $? -eq 0 ]]; then
            green "Bowtie2索引构建完成"
        else
            red "Bowtie2索引构建失败,请检查日志: $outdir/logs/bowtie2_build.log"
            exit 1
        fi
    else
        yellow "Bowtie2索引已存在,跳过构建"
    fi
    
    # 保存基因组信息
    cat > "$genome_dir/genome_info.txt" << EOF
# 基因组信息
GENOME_FILE=$genome_file
FRAGMENTS_FILE=$fragments_file
BOWTIE2_INDEX=$bt2_index_base
RESTRICTION_SITE=$restriction_site
EOF
    
    green "基因组预处理完成"
    
    # 返回重要文件路径
    echo "$fragments_file"
    echo "$bt2_index_base"
}

# 读取样本表函数
read_samplesheet() {
    local samplesheet="$1"
    local outdir="$2"
    
    blue "读取样本表: $samplesheet"
    
    # 验证样本表格式并转换
    local validated_samplesheet="$outdir/logs/samples.validated.csv"
    python3 d:\download\下载\hic-2.1.0\bin\check_samplesheet.py "$samplesheet" > "$validated_samplesheet"
    
    if [[ $? -ne 0 ]]; then
        red "样本表格式验证失败"
        exit 1
    fi
    
    # 解析样本信息
    local samples=()
    local fastq_files=()
    
    # 跳过标题行,读取样本数据
    tail -n +2 "$validated_samplesheet" | while IFS=, read -r sample fastq_1 fastq_2 single_end; do
        # 去除引号
        sample=$(echo "$sample" | tr -d '"')
        fastq_1=$(echo "$fastq_1" | tr -d '"')
        fastq_2=$(echo "$fastq_2" | tr -d '"')
        single_end=$(echo "$single_end" | tr -d '"')
        
        # 验证FASTQ文件存在
        if [[ ! -f "$fastq_1" ]]; then
            red "FASTQ文件不存在: $fastq_1"
            exit 1
        fi
        
        if [[ "$single_end" == "false" ]] && [[ ! -f "$fastq_2" ]]; then
            red "FASTQ文件不存在: $fastq_2"
            exit 1
        fi
        
        # 输出样本信息供主脚本使用
        echo "SAMPLE:$sample"
        echo "FASTQ1:$fastq_1"
        echo "FASTQ2:$fastq_2"
        echo "SINGLE_END:$single_end"
        echo "---"
    done
}

# 运行FASTQC
run_fastqc() {
    local fastq_file="$1"
    local outdir="$2"
    local sample_name="$3"
    
    blue "运行FASTQC: $sample_name"
    
    local fastqc_outdir="$outdir/fastqc"
    mkdir -p "$fastqc_outdir"
    
    # 检查是否已存在结果
    local fastqc_zip="$fastqc_outdir/$(basename "${fastq_file%.*}")_fastqc.zip"
    if [[ ! -f "$fastqc_zip" ]]; then
        fastqc \
            --threads 1 \
            --outdir "$fastqc_outdir" \
            --quiet \
            "$fastq_file" \
            > "$outdir/logs/fastqc_${sample_name}.log" 2>&1
        
        if [[ $? -eq 0 ]]; then
            green "FASTQC完成: $sample_name"
        else
            red "FASTQC失败: $sample_name"
            return 1
        fi
    else
        yellow "FASTQC结果已存在,跳过: $sample_name"
    fi
}

# 运行Bowtie2映射
run_bowtie2_mapping() {
    local fastq_1="$1"
    local fastq_2="$2"
    local bt2_index="$3"
    local outdir="$4"
    local sample_name="$5"
    local cpu="$6"
    local min_mapq="$7"
    
    blue "运行Bowtie2映射: $sample_name"
    
    local mapping_dir="$outdir/hicpro/mapping"
    mkdir -p "$mapping_dir"
    
    local output_bam="$mapping_dir/${sample_name}.bam"
    local output_log="$mapping_dir/${sample_name}.bowtie2.log"
    
    if [[ ! -f "$output_bam" ]]; then
        # 构建bowtie2命令
        local bowtie2_cmd="bowtie2"
        bowtie2_cmd+=" --very-sensitive"  # 使用高敏感度模式
        bowtie2_cmd+=" -L 30"            # 种子长度
        bowtie2_cmd+=" --score-min L,-0.6,-0.2"  # 评分阈值
        bowtie2_cmd+=" --end-to-end"     # 端到端对齐
        bowtie2_cmd+=" --reorder"        # 保持读取顺序
        bowtie2_cmd+=" --threads $cpu"
        bowtie2_cmd+=" -x $bt2_index"
        
        if [[ -n "$fastq_2" ]] && [[ "$fastq_2" != "null" ]]; then
            # 双端测序
            bowtie2_cmd+=" -1 $fastq_1 -2 $fastq_2"
        else
            # 单端测序
            bowtie2_cmd+=" -U $fastq_1"
        fi
        
        # 执行bowtie2映射并转换为BAM格式
        blue "执行Bowtie2映射..."
        $bowtie2_cmd \
            2> "$output_log" \
            | samtools view -bS -q "$min_mapq" - \
            | samtools sort -@ "$cpu" -o "$output_bam" - \
            > "$outdir/logs/bowtie2_${sample_name}.log" 2>&1
        
        if [[ $? -eq 0 ]]; then
            # 创建BAM索引
            samtools index "$output_bam"
            green "Bowtie2映射完成: $sample_name"
        else
            red "Bowtie2映射失败: $sample_name"
            return 1
        fi
    else
        yellow "Bowtie2映射结果已存在,跳过: $sample_name"
    fi
    
    echo "$output_bam"
}

# 获取有效互作
get_valid_interactions()
