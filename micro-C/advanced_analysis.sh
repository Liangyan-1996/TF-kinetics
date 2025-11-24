#!/bin/bash

# Hi-C高级分析功能模块
# 包含compartments和TADs分析功能

# 运行compartments分析
run_compartments_analysis() {
    local cooler_file="$1"
    local genome_fasta="$2"
    local chr_size_file="$3"
    local outdir="$4"
    local sample_name="$5"
    local resolution="$6"
    
    blue "运行compartments分析: $sample_name (分辨率: $resolution)"
    
    local compartments_dir="$outdir/compartments"
    mkdir -p "$compartments_dir"
    
    local compartments_file="$compartments_dir/${sample_name}_${resolution}_compartments.bedgraph"
    
    if [[ ! -f "$compartments_file" ]]; then
        if command -v cooltools >/dev/null 2>&1; then
            blue "使用cooltools进行compartments分析..."
            
            # 创建基因组bin文件
            local genome_bins="$compartments_dir/genome_bins_${resolution}.txt"
            cooltools genome binnify --all-names "$chr_size_file" "$resolution" > "$genome_bins" 2>"$outdir/logs/compartments_bins_${sample_name}.log"
            
            # 计算GC含量
            local gc_file="$compartments_dir/genome_gc_${resolution}.txt"
            cooltools genome gc "$genome_bins" "$genome_fasta" > "$gc_file" 2>"$outdir/logs/compartments_gc_${sample_name}.log"
            
            # 运行eigs-cis分析
            cooltools eigs-cis -o "${compartments_dir}/${sample_name}_${resolution}_compartments" "$cooler_file" 2>"$outdir/logs/compartments_${sample_name}.log"
            
            if [[ $? -eq 0 ]]; then
                green "compartments分析完成: $sample_name ($resolution)"
                
                # 转换为bigwig格式（如果可能）
                if command -v bedGraphToBigWig >/dev/null 2>&1; then
                    local bw_file="$compartments_dir/${sample_name}_${resolution}_compartments.bw"
                    bedGraphToBigWig "$compartments_file" "$chr_size_file" "$bw_file" 2>/dev/null || true
                fi
            else
                red "compartments分析失败: $sample_name ($resolution)"
                return 1
            fi
        else
            yellow "cooltools未安装,创建占位符文件"
            echo "# Compartments analysis placeholder" > "$compartments_file"
            echo "# cooltools not available" >> "$compartments_file"
        fi
    else
        yellow "compartments文件已存在,跳过: $sample_name ($resolution)"
    fi
    
    echo "$compartments_file"
}

# 运行TADs分析 - insulation方法
run_tads_insulation() {
    local cooler_file="$1"
    local outdir="$2"
    local sample_name="$3"
    local resolution="$4"
    local window_size="${5:-100000}"  # 默认窗口大小100kb
    
    blue "运行TADs分析 (insulation方法): $sample_name (分辨率: $resolution)"
    
    local tads_dir="$outdir/tads"
    mkdir -p "$tads_dir"
    
    local insulation_file="$tads_dir/${sample_name}_${resolution}_insulation.tsv"
    local boundaries_file="$tads_dir/${sample_name}_${resolution}_boundaries.bed"
    
    if [[ ! -f "$insulation_file" ]]; then
        if command -v cooltools >/dev/null 2>&1; then
            blue "使用cooltools insulation分析..."
            
            # 运行insulation分析
            cooltools insulation \
                --window "$window_size" \
                --min-dist-bad-bin 2 \
                "$cooler_file" \
                > "$insulation_file" 2>"$outdir/logs/tads_insulation_${sample_name}.log"
            
            if [[ $? -eq 0 ]]; then
                # 提取TAD边界
                awk 'NR>1 && $6>0.5 {print $1"\t"$2"\t"$3"\tTAD_boundary\t"$6}' "$insulation_file" > "$boundaries_file"
                
                green "TADs insulation分析完成: $sample_name ($resolution)"
            else
                red "TADs insulation分析失败: $sample_name ($resolution)"
                return 1
            fi
        else
            yellow "cooltools未安装,创建占位符文件"
            echo "# TADs insulation analysis placeholder" > "$insulation_file"
            echo "# cooltools not available" >> "$insulation_file"
            touch "$boundaries_file"
        fi
    else
        yellow "TADs insulation文件已存在,跳过: $sample_name ($resolution)"
    fi
    
    echo "$insulation_file"
    echo "$boundaries_file"
}

# 运行TADs分析 - HiCExplorer方法
run_tads_hicexplorer() {
    local cooler_file="$1"
    local outdir="$2"
    local sample_name="$3"
    local resolution="$4"
    local min_depth="${5:-0.01}"
    local max_depth="${6:-0.1}"
    
    blue "运行TADs分析 (HiCExplorer方法): $sample_name (分辨率: $resolution)"
    
    local tads_dir="$outdir/tads"
    mkdir -p "$tads_dir"
    
    local tads_file="$tads_dir/${sample_name}_${resolution}_tads.bed"
    local boundaries_file="$tads_dir/${sample_name}_${resolution}_tad_boundaries.bed"
    
    if [[ ! -f "$tads_file" ]]; then
        if command -v hicFindTADs >/dev/null 2>&1; then
            blue "使用hicFindTADs分析..."
            
            # 转换为h5格式（如果需要）
            local h5_file="$tads_dir/${sample_name}_${resolution}.h5"
            if [[ "$cooler_file" =~ \.cool$ ]]; then
                if command -v cooler2hic5 >/dev/null 2>&1; then
                    cooler2hic5 "$cooler_file" "$h5_file" 2>/dev/null || cp "$cooler_file" "$h5_file"
                else
                    cp "$cooler_file" "$h5_file"
                fi
            else
                cp "$cooler_file" "$h5_file"
            fi
            
            # 运行TAD检测
            hicFindTADs \
                --matrix "$h5_file" \
                --outPrefix "$tads_dir/${sample_name}_${resolution}" \
                --minDepth "$min_depth" \
                --maxDepth "$max_depth" \
                --numberOfProcessors 1 \
                > "$outdir/logs/tads_hicexplorer_${sample_name}.log" 2>&1
            
            if [[ $? -eq 0 ]]; then
                # 整理输出文件
                find "$tads_dir" -name "*tad_score*" -type f | head -n1 | xargs -I {} cp {} "$tads_file"
                find "$tads_dir" -name "*boundaries*" -type f | head -n1 | xargs -I {} cp {} "$boundaries_file"
                
                green "TADs HiCExplorer分析完成: $sample_name ($resolution)"
            else
                red "TADs HiCExplorer分析失败: $sample_name ($resolution)"
                return 1
            fi
        else
            yellow "hicFindTADs未安装,创建占位符文件"
            echo "# TADs HiCExplorer analysis placeholder" > "$tads_file"
            echo "# hicFindTADs not available" >> "$tads_file"
            touch "$boundaries_file"
        fi
    else
        yellow "TADs HiCExplorer文件已存在,跳过: $sample_name ($resolution)"
    fi
    
    echo "$tads_file"
    echo "$boundaries_file"
}

# 整合TADs分析结果
merge_tads_results() {
    local outdir="$1"
    local method="$2"  # insulation 或 hicexplorer
    
    blue "整合TADs分析结果: $method"
    
    local tads_dir="$outdir/tads"
    local merged_file="$tads_dir/merged_tads_${method}.bed"
    
    # 收集所有TAD边界文件
    local boundary_files=()
    if [[ "$method" == "insulation" ]]; then
        mapfile -t boundary_files < <(find "$tads_dir" -name "*_insulation.tsv" -type f)
    else
        mapfile -t boundary_files < <(find "$tads_dir" -name "*_tad_boundaries.bed" -type f)
    fi
    
    if [[ ${#boundary_files[@]} -gt 0 ]]; then
        # 合并边界文件
        {
            echo -e "#chrom\tstart\tend\tsample\tscore"
            for boundary_file in "${boundary_files[@]}"; do
                local sample_name=$(basename "$boundary_file" | sed "s/_${method}.*//")
                if [[ -s "$boundary_file" ]]; then
                    awk -v sample="$
