#!/bin/bash

# Hi-C Pipeline 使用示例和测试脚本

# 设置基础参数
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOME_FILE="$SCRIPT_DIR/test_data/genome.fa"
SAMPLES_FILE="$SCRIPT_DIR/test_data/samples.csv"
OUTDIR="$SCRIPT_DIR/test_results"

# 创建测试数据目录
mkdir -p "$SCRIPT_DIR/test_data"
mkdir -p "$OUTDIR"

# 创建测试基因组文件（简化版本）
create_test_genome() {
    echo "创建测试基因组文件..."
    cat > "$GENOME_FILE" << 'EOF'
>chr1
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
>chr2
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
    
    # 创建染色体大小文件
    cat > "$SCRIPT_DIR/test_data/chrom_sizes.txt" << 'EOF'
chr1	186
nchr2	186
EOF
}

# 创建测试样本表
create_test_samples() {
    echo "创建测试样本表..."
    cat > "$SAMPLES_FILE" << 'EOF'
sample,fastq_1,fastq_2
sample1,$SCRIPT_DIR/test_data/sample1_R1.fastq,$SCRIPT_DIR/test_data/sample1_R2.fastq
sample2,$SCRIPT_DIR/test_data/sample2_R1.fastq,$SCRIPT_DIR/test_data/sample2_R2.fastq
EOF
}

# 创建测试FASTQ文件（简化版本）
create_test_fastq() {
    echo "创建测试FASTQ文件..."
    
    # Sample 1 - Read 1
    cat > "$SCRIPT_DIR/test_data/sample1_R1.fastq" << 'EOF'
@read1_sample1
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_sample1
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    
    # Sample 1 - Read 2
    cat > "$SCRIPT_DIR/test_data/sample1_R2.fastq" << 'EOF'
@read1_sample1
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_sample1
CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    
    # Sample 2 - Read 1
    cat > "$SCRIPT_DIR/test_data/sample2_R1.fastq" << 'EOF'
@read1_sample2
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_sample2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
    
    # Sample 2 - Read 2
    cat > "$SCRIPT_DIR/test_data/sample2_R2.fastq" << 'EOF'
@read1_sample2
TACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2_sample2
CGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF
}

# 创建酶切位点文件
create_resfrag_file() {
    echo "创建酶切片段文件..."
    
    # 使用digest_genome.py创建酶切片段文件
    if [[ -f "$SCRIPT_DIR/bin/digest_genome.py" ]]; then
        python3 "$SCRIPT_DIR/bin/digest_genome.py" \
            "$GENOME_FILE" \
            "A^AGCTT" \
            "$OUTDIR/genome/resfrag.bed" \
            2>"$OUTDIR/logs/digest_genome.log" || true
    else
        # 创建占位符酶切片段文件
        mkdir -p "$OUTDIR/genome"
        cat > "$OUTDIR/genome/resfrag.bed" << 'EOF'
chr1	0	50	HindIII_fragment_1	0	+
chr1	50	100	HindIII_fragment_2	0	+
chr1	100	150	HindIII_fragment_3	0	+
chr1	150	186	HindIII_fragment_4	0	+
chr2	0	50	HindIII_fragment_5	0	+
chr2	50	100	HindIII_fragment_6	0	+
chr2	100	150	HindIII_fragment_7	0	+
chr2	150	186	HindIII_fragment_8	0	+
EOF
    fi
}

# 运行依赖检查
run_dependency_check() {
    echo "运行依赖检查..."
    bash "$SCRIPT_DIR/check_dependencies.sh" --output "$OUTDIR/logs/dependency_check.log"
}

# 运行完整pipeline测试
run_pipeline_test() {
    echo "运行完整pipeline测试..."
    
    local main_script="$SCRIPT_DIR/hic_pipeline.sh"
    
    if [[ -f "$main_script" ]]; then
        bash "$main_script" \
            -g "$GENOME_FILE" \
            -s "$SAMPLES_FILE" \
            -r "A^AGCTT" \
            -l "AGCT" \
            -o "$OUTDIR" \
            --resolutions "100000,50000,25000" \
            --cpu 2 \
            --memory "4G" \
            --run_advanced \
            --run_compartments \
            --run_tads \
            --tads_method "insulation" \
            2>&1 | tee "$OUTDIR/logs/pipeline_test.log"
    else
        echo "错误: 主脚本文件不存在: $main_script"
        return 1
    fi
}

# 验证输出结果
validate_results() {
    echo "验证输出结果..."
    
    local validation_passed=true
    
    # 检查主要输出目录
    local expected_dirs=(
        "$OUTDIR/fastqc"
        "$OUTDIR/mapping"
        "$OUTDIR/hicpro"
        "$OUTDIR/analysis"
        "$OUTDIR/plots"
        "$OUTDIR/reports"
        "$OUTDIR/compartments"
        "$OUTDIR/tads"
    )
    
    for dir in "${expected_dirs[@]}"; do
        if [[ -d "$dir" ]]; then
            echo "✓ 目录存在: $(basename "$dir")"
        else
            echo "✗ 目录缺失: $(basename "$dir")"
            validation_passed=false
        fi
    done
    
    # 检查关键文件
    local expected_files=(
        "$OUTDIR/genome/resfrag.bed"
        "$OUTDIR/genome/bowtie2_index.1.bt2"
        "$OUTDIR/analysis/merged_samples.cool"
        "$OUTDIR/reports/final_report.html"
    )
    
    for file in "${expected_files[@]}"; do
        if [[ -f "$file" ]]; then
            echo "✓ 文件存在: $(basename "$file")"
        else
            echo "✗ 文件缺失: $(basename "$file")"
            validation_passed=false
        fi
    done
    
    if [[ "$validation_passed" == "true" ]]; then
        echo "✓ 所有验证通过！"
    else
        echo "✗ 验证失败，请检查日志文件"
    fi
    
    return $([ "$validation_passed" == "true" ] && echo 0 || echo 1)
}

# 生成测试报告
generate_test_report() {
    echo "生成测试报告..."
    
    local report_file="$OUTDIR/test_report.html"
    
    cat > "$report_file" << 'EOF'
<!DOCTYPE html>
<html>
<head>
    <title>Hi-C Pipeline 测试报告</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        .header { background-color: #4CAF50; color: white; padding: 10px; }
        .section { margin: 20px 0; }
        .success { color: green; }
        .error { color: red; }
        .warning { color: orange; }
        table { border-collapse: collapse; width: 100%; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Hi-C Pipeline 测试报告</h1>
        <p>生成时间: $(date)</p>
    </div>
    
    <div class="section">
        <h2>测试概述</h2>
        <p>本测试验证了Hi-C数据分析pipeline的基本功能，包括：</p>
        <ul>
            <li>基因组预处理（酶切片段生成）</li>
            <li>FASTQC质量控制</li>
            <li>Bowtie2序列比对</li>
            <li>Hi-C互作分析</li>
            <li>接触矩阵构建</li>
            <li>高级分析（compartments和TADs）</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>输出文件结构</h2>
        <table>
            <tr><th>目录</th><th>内容</th></tr>
            <tr><td>fastqc/</td><td>FASTQC质量控制报告</td></tr>
            <tr><td>mapping/</td><td>Bowtie2比对结果</td></tr>
            <tr><td>hicpro/</td><td>Hi-C互作分析结果</td></tr>
            <tr><td>analysis/</td><td>接触矩阵和分析结果</td></tr>
            <tr><td>plots/</td><td>可视化图表</td></tr>
            <tr><td>compartments/</td><td>Compartments分析结果</td></tr>
            <tr><td>tads/</td><td>TADs分析结果</td></tr>
            <tr><td>reports/</td><td>最终分析报告</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>使用说明</h2>
        <p>要使用此pipeline分析您的Hi-C数据，请运行：</p>
        <pre><code>bash hic_pipeline.sh -g your_genome.fa -s your_samples.csv -r A^AGCTT -l AGCT -o output_dir</code></pre>
        
        <p>要运行高级分析（compartments和TADs），添加：</p>
        <pre><code>--run_advanced --run_compartments --run_tads --tads_method insulation</code></pre>
    </div>
    
    <div class="section">
        <h2>故障排除</h2>
        <p>如果遇到问题，请检查：</p>
        <ul>
            <li>依赖软件是否安装完整（运行 check_dependencies.sh）</li>
            <li>输入文件格式是否正确</li>
            <li>日志文件中的错误信息</li>
            <li>系统资源是否充足</li>
        </ul>
    </div>
</body>
</html>
EOF
    
    echo "测试报告已生成: $report_file"
}

# 主函数
main() {
    echo "=== Hi-C Pipeline 测试脚本 ==="
    echo "脚本目录: $SCRIPT_DIR"
    echo "输出目录: $OUTDIR"
    echo
    
    # 创建测试数据
    create_test_genome
    create_test_samples
    create_test_fastq
    create_resfrag_file
    
    # 运行依赖检查
    run_dependency_check
    
    # 运行pipeline测试
    if run_pipeline_test; then
        echo "Pipeline测试运行完成"
    else
        echo "Pipeline测试运行失败"
    fi
    
    # 验证结果
    validate_results
    
    # 生成测试报告
    generate_test_report
    
    echo
    echo "=== 测试完成 ==="
    echo "查看详细结果: $OUTDIR"
    echo "查看测试报告: $report_file"
}

# 如果直接运行此脚本
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi