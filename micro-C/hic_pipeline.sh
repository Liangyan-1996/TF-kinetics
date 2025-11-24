#!/bin/bash

# Hi-C数据分析Pipeline (Shell脚本版本)
# 基于nf-core/hic pipeline转换而来, 用于无容器环境下的Hi-C数据分析

set -euo pipefail

# 默认参数
GENOME=""
SAMPLESHEET=""
RESTRICTION_SITE=""
LIGATION_SITE=""
OUTDIR="./results"
RESOLUTIONS="1000000,500000,250000,100000,50000,25000,10000,5000,1000"
MIN_MAPQ=20
MAX_INSERT_SIZE=2000
MIN_INSERT_SIZE=50
MIN_FRAG_SIZE=50
MAX_FRAG_SIZE=100000
MIN_CIS_DIST=1000
NCPU=4
MEMORY="8G"
RUN_ADVANCED=false
RUN_COMPARTMENTS=false
RUN_TADS=false
TADS_METHOD="insulation"

# 颜色输出函数
red() { echo -e "\033[31m$1\033[0m"; }
green() { echo -e "\033[32m$1\033[0m"; }
yellow() { echo -e "\033[33m$1\033[0m"; }
blue() { echo -e "\033[34m$1\033[0m"; }

# 帮助信息
usage() {
    cat << EOF
用法: $0 [选项]

Hi-C数据分析Pipeline

必需参数:
    -g, --genome <文件>          参考基因组FASTA文件
    -s, --samplesheet <文件>     样本表文件(CSV格式, 包含sample,fastq_1,fastq_2列)
    -r, --restriction_site <位点>  限制酶切位点(如: A^AGCTT)
    -l, --ligation_site <位点>   连接位点(如: AGCT)

可选参数:
    -o, --outdir <目录>          输出目录(默认: ./results)
    --resolutions <分辨率>       分析分辨率, 逗号分隔(默认: 1000000,500000,250000,100000,50000,25000,10000,5000,1000)
    --min_mapq <数值>             最小映射质量(默认: 20)
    --max_insert_size <数值>       最大插入片段大小(默认: 2000)
    --min_insert_size <数值>       最小插入片段大小(默认: 50)
    --min_frag_size <数值>        最小片段大小(默认: 50)
    --max_frag_size <数值>        最大片段大小(默认: 100000)
    --min_cis_dist <数值>         最小顺式距离(默认: 1000)
    --cpu <数值>                  CPU线程数(默认: 4)
    --memory <内存>               内存使用(默认: 8G)
    --run_advanced               运行高级分析(compartments和TADs)
    --run_compartments           运行compartments分析
    --run_tads                   运行TADs分析
    --tads_method <方法>         TADs分析方法(insulation或hicexplorer, 默认: insulation)
    -h, --help                    显示帮助信息

示例:
    $0 -g genome.fa -s samples.csv -r A^AGCTT -l AGCT -o output_dir

EOF
}

# 参数解析
parse_args() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -g|--genome)
                GENOME="$2"
                shift 2
                ;;
            -s|--samplesheet)
                SAMPLESHEET="$2"
                shift 2
                ;;
            -r|--restriction_site)
                RESTRICTION_SITE="$2"
                shift 2
                ;;
            -l|--ligation_site)
                LIGATION_SITE="$2"
                shift 2
                ;;
            -o|--outdir)
                OUTDIR="$2"
                shift 2
                ;;
            --resolutions)
                RESOLUTIONS="$2"
                shift 2
                ;;
            --min_mapq)
                MIN_MAPQ="$2"
                shift 2
                ;;
            --max_insert_size)
                MAX_INSERT_SIZE="$2"
                shift 2
                ;;
            --min_insert_size)
                MIN_INSERT_SIZE="$2"
                shift 2
                ;;
            --min_frag_size)
                MIN_FRAG_SIZE="$2"
                shift 2
                ;;
            --max_frag_size)
                MAX_FRAG_SIZE="$2"
                shift 2
                ;;
            --min_cis_dist)
                MIN_CIS_DIST="$2"
                shift 2
                ;;
            --cpu)
                NCPU="$2"
                shift 2
                ;;
            --memory)
                MEMORY="$2"
                shift 2
                ;;
            --run_advanced)
                RUN_ADVANCED=true
                shift
                ;;
            --run_compartments)
                RUN_COMPARTMENTS=true
                shift
                ;;
            --run_tads)
                RUN_TADS=true
                shift
                ;;
            --tads_method)
                TADS_METHOD="$2"
                shift 2
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                red "未知参数: $1"
                usage
                exit 1
                ;;
        esac
    done
}

# 检查依赖
check_dependencies() {
    blue "检查依赖软件..."
    
    local missing_deps=()
    
    # 检查基础工具
    command -v python3 >/dev/null 2>&1 || missing_deps+=("python3")
    command -v samtools >/dev/null 2>&1 || missing_deps+=("samtools")
    command -v bowtie2 >/dev/null 2>&1 || missing_deps+=("bowtie2")
    
    # 检查Python包
    python3 -c "import pysam" >/dev/null 2>&1 || missing_deps+=("python3-pysam")
    python3 -c "import bx" >/dev/null 2>&1 || missing_deps+=("python3-bx-python")
    
    if [[ ${#missing_deps[@]} -gt 0 ]]; then
        red "缺少依赖: ${missing_deps[*]}"
        yellow "请安装缺失的依赖后再运行此脚本"
        exit 1
    fi
    
    green "所有依赖检查通过"
}

# 检查输入文件
validate_inputs() {
    blue "验证输入文件..."
    
    # 检查必需文件
    if [[ -z "$GENOME" ]]; then
        red "错误: 必须指定参考基因组文件 (-g)"
        usage
        exit 1
    fi
    
    if [[ -z "$SAMPLESHEET" ]]; then
        red "错误: 必须指定样本表文件 (-s)"
        usage
        exit 1
    fi
    
    if [[ -z "$RESTRICTION_SITE" ]]; then
        red "错误: 必须指定限制酶切位点 (-r)"
        usage
        exit 1
    fi
    
    if [[ -z "$LIGATION_SITE" ]]; then
        red "错误: 必须指定连接位点 (-l)"
        usage
        exit 1
    fi
    
    # 检查文件是否存在
    if [[ ! -f "$GENOME" ]]; then
        red "错误: 参考基因组文件不存在: $GENOME"
        exit 1
    fi
    
    if [[ ! -f "$SAMPLESHEET" ]]; then
        red "错误: 样本表文件不存在: $SAMPLESHEET"
        exit 1
    fi
    
    # 验证样本表格式
    if ! python3 d:\download\下载\hic-2.1.0\bin\check_samplesheet.py "$SAMPLESHEET"; then
        red "错误: 样本表格式验证失败"
        exit 1
    fi
    
    green "输入文件验证通过"
}

# 创建输出目录结构
create_output_structure() {
    blue "创建输出目录结构..."
    
    mkdir -p "$OUTDIR"/{logs,genome,fastqc,hicpro,cooler,plots,multiqc,reports}
    mkdir -p "$OUTDIR"/hicpro/{mapping,valid_pairs,pairs,raw_maps,iced_maps}
    
    green "输出目录结构创建完成"
}

# 主函数
main() {
    echo "========================================"
    blue "Hi-C数据分析Pipeline (Shell脚本版本)"
    blue "基于nf-core/hic pipeline转换"
    echo "========================================"
    
    # 解析参数
    parse_args "$@"
    
    # 检查依赖
    check_dependencies
    
    # 验证输入
    validate_inputs
    
    # 创建输出目录
    create_output_structure
    
    echo ""
    green "Pipeline初始化完成"
    yellow "参考基因组: $GENOME"
    yellow "样本表: $SAMPLESHEET"
    yellow "限制酶切位点: $RESTRICTION_SITE"
    yellow "连接位点: $LIGATION_SITE"
    yellow "输出目录: $OUTDIR"
    yellow "分辨率: $RESOLUTIONS"
    echo ""
    
    # 这里将继续添加具体的处理步骤
    blue "开始执行分析步骤..."
    
    # 记录开始时间
    local start_time=$(date +%s)
    
    # 加载功能模块
    source "$(dirname "$0")/genome_functions.sh"
    source "$(dirname "$0")/analysis_functions.sh"
    source "$(dirname "$0")/advanced_analysis.sh"
    
    # 步骤1: 基因组预处理
    blue "步骤1: 基因组预处理"
    local genome_info=$(prepare_genome "$GENOME" "$RESTRICTION_SITE" "$OUTDIR" "$NCPU")
    local fragments_file=$(echo "$genome_info" | head -n1)
    local bt2_index=$(echo "$genome_info" | tail -n1)
    
    # 步骤2: FASTQC质量控制
    blue "步骤2: FASTQC质量控制"
    while IFS= read -r line; do
        if [[ "$line" =~ ^SAMPLE:(.*) ]]; then
            sample_name="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^FASTQ1:(.*) ]]; then
            fastq_1="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^FASTQ2:(.*) ]]; then
            fastq_2="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^SINGLE_END:(.*) ]]; then
            single_end="${BASH_REMATCH[1]}"
        elif [[ "$line" == "---" ]]; then
            # 处理完一个样本
            run_fastqc "$fastq_1" "$OUTDIR" "$sample_name"
            if [[ "$single_end" == "false" ]] && [[ "$fastq_2" != "null" ]]; then
                run_fastqc "$fastq_2" "$OUTDIR" "${sample_name}_R2"
            fi
        fi
    done < <(read_samplesheet "$SAMPLESHEET" "$OUTDIR")
    
    # 步骤3: HICPRO核心处理流程
    blue "步骤3: HICPRO核心处理流程"
    declare -a cooler_files=()
    declare -a valid_pairs_files=()
    
    while IFS= read -r line; do
        if [[ "$line" =~ ^SAMPLE:(.*) ]]; then
            sample_name="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^FASTQ1:(.*) ]]; then
            fastq_1="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^FASTQ2:(.*) ]]; then
            fastq_2="${BASH_REMATCH[1]}"
        elif [[ "$line" =~ ^SINGLE_END:(.*) ]]; then
            single_end="${BASH_REMATCH[1]}"
        elif [[ "$line" == "---" ]]; then
            # 处理完一个样本
            
            # Bowtie2映射
            local bam_file=$(run_bowtie2_mapping "$fastq_1" "$fastq_2" "$bt2_index" "$OUTDIR" "$sample_name" "$NCPU" "$MIN_MAPQ")
            
            # 获取有效互作
            local valid_pairs_file=$(get_valid_interactions "$bam_file" "$fragments_file" "$OUTDIR" "$sample_name" "$MIN_INSERT_SIZE" "$MAX_INSERT_SIZE" "$MIN_FRAG_SIZE" "$MAX_FRAG_SIZE" "$MIN_CIS_DIST")
            valid_pairs_files+=("$valid_pairs_file")
            
            # 对每个分辨率构建矩阵
            IFS=',' read -ra RES_ARRAY <<< "$RESOLUTIONS"
            for resolution in "${RES_ARRAY[@]}"; do
                # 构建接触矩阵
                local matrix_file=$(build_contact_matrix "$valid_pairs_file" "$fragments_file" "$resolution" "$OUTDIR" "$sample_name")
                
                # ICE标准化
                local iced_matrix=$(run_ice_normalization "$matrix_file" "$OUTDIR" "$sample_name" "$resolution")
                
                # 转换为COOLER格式
                local cooler_file=$(convert_to_cooler "$iced_matrix" "$fragments_file" "$resolution" "$OUTDIR" "$sample_name")
                cooler_files+=("$cooler_file")
            done
        fi
    done < <(read_samplesheet "$SAMPLESHEET" "$OUTDIR")
    
    # 步骤4: 合并COOLER文件
    blue "步骤4: 合并COOLER文件"
    IFS=',' read -ra RES_ARRAY <<< "$RESOLUTIONS"
    for resolution in "${RES_ARRAY[@]}"; do
        # 收集该分辨率的所有cooler文件
        local res_cooler_files=()
        for cooler_file in "${cooler_files[@]}"; do
            if [[ "$cooler_file" =~ ${resolution}\.cool$ ]]; then
                res_cooler_files+=("$cooler_file")
            fi
        done
        
        if [[ ${#res_cooler_files[@]} -gt 0 ]]; then
            merge_cooler_files "${res_cooler_files[@]}"
        fi
    done
    
    # 步骤5: 可视化图表生成
    blue "步骤5: 可视化图表生成"
    for valid_pairs_file in "${valid_pairs_files[@]}"; do
        local sample_name=$(basename "$valid_pairs_file" .validPairs)
        plot_distance_vs_counts "$valid_pairs_file" "$OUTDIR" "$sample_name"
    done
    
    # 步骤6: MultiQC报告
    blue "步骤6: MultiQC报告生成"
    run_multiqc "$OUTDIR"
    
    # 记录结束时间
    local end_time=$(date +%s)
    
    # 步骤8: 高级分析(compartments和TADs)
    if [[ "${RUN_ADVANCED:-false}" == "true" ]]; then
        blue "步骤8: 高级分析(compartments和TADs)"
        
        # 对每个样本运行高级分析
        while IFS= read -r line; do
            if [[ "$line" =~ ^SAMPLE:(.*) ]]; then
                sample_name="${BASH_REMATCH[1]}"
            elif [[ "$line" =~ ^FASTQ1:(.*) ]]; then
                fastq_1="${BASH_REMATCH[1]}"
            elif [[ "$line" =~ ^FASTQ2:(.*) ]]; then
                fastq_2="${BASH_REMATCH[1]}"
            elif [[ "$line" =~ ^SINGLE_END:(.*) ]]; then
                single_end="${BASH_REMATCH[1]}"
            elif [[ "$line" == "---" ]]; then
                # 对每个分辨率运行分析
                IFS=',' read -ra RES_ARRAY <<< "$RESOLUTIONS"
                for resolution in "${RES_ARRAY[@]}"; do
                    local cooler_file="$OUTDIR/analysis/${sample_name}_${resolution}.cool"
                    
                    if [[ -f "$cooler_file" ]]; then
                        run_advanced_analysis \
                            "$cooler_file" \
                            "$GENOME" \
                            "$OUTDIR" \
                            "$sample_name" \
                            "$resolution" \
                            "${RUN_COMPARTMENTS:-true}" \
                            "${RUN_TADS:-true}" \
                            "${TADS_METHOD:-insulation}"
                    else
                        yellow "cooler文件不存在, 跳过高级分析: $sample_name ($resolution)"
                    fi
                done
            fi
        done < <(read_samplesheet "$SAMPLESHEET" "$OUTDIR")
        
        # 整合TADs结果
        if [[ "${RUN_TADS:-true}" == "true" ]]; then
            merge_tads_results "$OUTDIR" "${TADS_METHOD:-insulation}"
        fi
    fi
    
    # 生成最终报告
    blue "步骤9: 生成最终分析报告"
    
    local end_time=$(date +%s)
    local total_time=$((end_time - start_time))
    
    generate_final_report "$OUTDIR" "$SAMPLESHEET" "$GENOME" "$RESFRAG" "$total_time"
    
    echo ""
    green "========================================"
    green "Hi-C数据分析Pipeline执行完成"
    green "总运行时间: $((end_time - start_time))秒"
    green "结果目录: $OUTDIR"
    green "========================================"
}

# 如果脚本被直接执行, 则运行主函数
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi