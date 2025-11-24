#!/bin/bash

# 项目概览脚本 - 显示所有创建的文件和功能

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "========================================"
echo "Hi-C Pipeline Shell版本 - 项目概览"
echo "========================================"
echo

echo "📁 主要脚本文件:"
echo "├── hic_pipeline.sh               - 主pipeline脚本"
echo "├── check_dependencies.sh         - 依赖检查脚本"
echo "├── genome_functions.sh           - 基因组处理功能"
echo "├── analysis_functions.sh         - 分析功能模块"
echo "├── advanced_analysis.sh          - 高级分析功能"
echo "├── run_test.sh                   - 测试脚本"
echo "└── README_SHELL_VERSION.md       - 使用文档"
echo

echo "🔧 功能模块概览:"
echo
echo "1. 基因组预处理模块 (genome_functions.sh):"
echo "   ✓ digest_genome()              - 基因组酶切片段生成"
echo "   ✓ build_bowtie2_index()        - Bowtie2索引构建"
echo "   ✓ run_fastqc()                 - FASTQC质量控制"
echo "   ✓ run_bowtie2_mapping()        - Bowtie2序列比对"
echo "   ✓ get_valid_interactions()     - 有效互作提取"
echo
echo "2. 分析功能模块 (analysis_functions.sh):"
echo "   ✓ build_contact_matrices()      - 接触矩阵构建"
echo "   ✓ run_ice_normalization()       - ICE标准化"
echo "   ✓ convert_to_cooler()           - COOLER格式转换"
echo "   ✓ merge_cooler_files()          - COOLER文件合并"
echo "   ✓ plot_distance_vs_counts()     - 距离vs计数图"
echo "   ✓ run_multiqc()                 - MultiQC报告生成"
echo "   ✓ generate_final_report()       - 最终报告生成"
echo
echo "3. 高级分析模块 (advanced_analysis.sh):"
echo "   ✓ run_compartments_analysis()   - Compartments分析"
echo "   ✓ run_tads_insulation()         - TADs分析(insulation)"
echo "   ✓ run_tads_hicexplorer()        - TADs分析(HiCExplorer)"
echo "   ✓ merge_tads_results()          - TADs结果整合"
echo "   ✓ plot_tads_visualization()     - TADs可视化"
echo "   ✓ run_advanced_analysis()       - 高级分析主函数"
echo
echo "📊 支持的酶切位点:"
echo "   • HindIII: A^AGCTT"
echo "   • MboI: GATC"
echo "   • DpnII: GATC"
echo "   • HinfI: G^ANTC"
echo "   • 支持自定义酶切位点"
echo
echo "🔍 分析分辨率:"
echo "   默认: 1Mb, 500kb, 250kb, 100kb, 50kb, 25kb, 10kb, 5kb, 1kb"
echo "   可自定义分辨率范围"
echo
echo "📈 输出文件类型:"
echo "   • BED/BEDGRAPH: 基因组坐标文件"
echo "   • COOLER: Hi-C接触矩阵格式"
echo "   • TSV: 表格数据"
echo "   • PNG: 可视化图表"
echo "   • HTML: 分析报告"
echo "   • BigWig: 基因组浏览器格式"
echo
echo "🚀 快速开始命令:"
echo "   # 基本分析"
echo "   bash hic_pipeline.sh -g genome.fa -s samples.csv -r A^AGCTT -l AGCT -o output"
echo
echo "   # 高级分析(包含compartments和TADs)"
echo "   bash hic_pipeline.sh -g genome.fa -s samples.csv -r A^AGCTT -l AGCT -o output \\"
echo "       --resolutions \"1000000,500000,250000\" --run_advanced --run_compartments --run_tads"
echo
echo "🧪 测试运行:"
echo "   bash run_test.sh"
echo
echo "📋 参数检查:"
echo "   bash check_dependencies.sh"
echo
echo "📖 详细文档:"
echo "   cat README_SHELL_VERSION.md"
echo
echo "========================================"
echo "状态: ✅ 所有功能模块已完成"
echo "基于: nf-core/hic v2.1.0"
echo "转换: Shell脚本版本"
echo "========================================"