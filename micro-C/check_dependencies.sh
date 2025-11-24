#!/bin/bash

# Hi-C Pipeline依赖检查和环境配置脚本

set -euo pipefail

# 颜色输出函数
red() { echo -e "\033[31m$1\033[0m"; }
green() { echo -e "\033[32m$1\033[0m"; }
yellow() { echo -e "\033[33m$1\033[0m"; }
blue() { echo -e "\033[34m$1\033[0m"; }

# 依赖列表
declare -a REQUIRED_COMMANDS=(
    "python3"
    "samtools"
    "bowtie2"
    "bowtie2-build"
    "fastqc"
)

declare -a PYTHON_PACKAGES=(
    "pysam"
    "bx-python"
    "numpy"
    "pandas"
    "matplotlib"
)

declare -a OPTIONAL_COMMANDS=(
    "cooler"
    "multiqc"
)

# 检查命令是否存在
check_command() {
    local cmd="$1"
    local optional="${2:-false}"
    
    if command -v "$cmd" >/dev/null 2>&1; then
        local version=$(get_version "$cmd")
        if [[ "$optional" == "true" ]]; then
            green "✓ $cmd (optional) - $version"
        else
            green "✓ $cmd - $version"
        fi
        return 0
    else
        if [[ "$optional" == "true" ]]; then
            yellow "⚠ $cmd (optional) - not installed"
        else
            red "✗ $cmd - not installed"
        fi
        return 1
    fi
}

# 获取软件版本
get_version() {
    local cmd="$1"
    local version=""
    
    case "$cmd" in
        "python3")
            version=$(python3 --version 2>&1 | head -n1)
            ;;
        "samtools")
            version=$(samtools --version 2>&1 | head -n1)
            ;;
        "bowtie2"|"bowtie2-build")
            version=$(bowtie2 --version 2>&1 | head -n1 | cut -d' ' -f3)
            ;;
        "fastqc")
            version=$(fastqc --version 2>&1 | head -n1)
            ;;
        "cooler")
            version=$(cooler --version 2>&1 | head -n1)
            ;;
        "multiqc")
            version=$(multiqc --version 2>&1 | head -n1)
            ;;
        *)
            version="version unknown"
            ;;
    esac
    
    echo "$version"
}

# 检查Python包
check_python_package() {
    local package="$1"
    local optional="${2:-false}"
    
    if python3 -c "import $package" >/dev/null 2>&1; then
        local version=$(python3 -c "import $package; print(getattr($package, '__version__', 'version unknown'))" 2>/dev/null || echo "version unknown")
        if [[ "$optional" == "true" ]]; then
            green "✓ Python package: $package (optional) - $version"
        else
            green "✓ Python package: $package - $version"
        fi
        return 0
    else
        if [[ "$optional" == "true" ]]; then
            yellow "⚠ Python package: $package (optional) - not installed"
        else
            red "✗ Python package: $package - not installed"
        fi
        return 1
    fi
}

# 检查系统资源
check_system_resources() {
    blue "\nChecking system resources..."
    
    # CPU信息
    if [[ -f /proc/cpuinfo ]]; then
        local cpu_count=$(grep -c "^processor" /proc/cpuinfo)
        local cpu_model=$(grep "model name" /proc/cpuinfo | head -n1 | cut -d':' -f2 | sed 's/^[[:space:]]*//')
        green "✓ CPU: $cpu_count cores - $cpu_model"
    else
        yellow "⚠ Unable to get CPU info"
    fi
    
    # 内存信息
    if command -v free >/dev/null 2>&1; then
        local total_mem=$(free -h | grep "^Mem:" | awk '{print $2}')
        local available_mem=$(free -h | grep "^Mem:" | awk '{print $7}')
        green "✓ Memory: total $total_mem, available $available_mem"
    else
        yellow "⚠ Unable to get memory info"
    fi
    
    # 磁盘空间
    local current_dir=$(pwd)
    if command -v df >/dev/null 2>&1; then
        local disk_info=$(df -h "$current_dir" | tail -n1)
        local disk_size=$(echo "$disk_info" | awk '{print $2}')
        local disk_used=$(echo "$disk_info" | awk '{print $3}')
        local disk_avail=$(echo "$disk_info" | awk '{print $4}')
        green "✓ Disk space: total $disk_size, used $disk_used, available $disk_avail"
    else
        yellow "⚠ Unable to get disk space info"
    fi
}

# 安装建议
show_installation_suggestions() {
    local missing_deps=()
    local missing_packages=()
    
    # 收集缺失的依赖
    for cmd in "${REQUIRED_COMMANDS[@]}"; do
        if ! check_command "$cmd" >/dev/null 2>&1; then
            missing_deps+=("$cmd")
        fi
    done
    
    # 收集缺失的Python包
    for pkg in "${PYTHON_PACKAGES[@]}"; do
        if ! check_python_package "$pkg" >/dev/null 2>&1; then
            missing_packages+=("$pkg")
        fi
    done
    
    if [[ ${#missing_deps[@]} -gt 0 ]] || [[ ${#missing_packages[@]} -gt 0 ]]; then
        blue "\nInstallation suggestions:"
        
        if [[ ${#missing_deps[@]} -gt 0 ]]; then
            echo "Missing command-line tools:"
            echo "# Ubuntu/Debian:"
            echo "sudo apt-get update"
            echo "sudo apt-get install -y ${missing_deps[*]}"
            echo ""
            echo "# CentOS/RHEL:"
            echo "sudo yum install -y ${missing_deps[*]}"
            echo ""
            echo "# macOS (using Homebrew):"
            echo "brew install ${missing_deps[*]}"
            echo ""
        fi
        
        if [[ ${#missing_packages[@]} -gt 0 ]]; then
            echo "Missing Python packages:"
            echo "pip3 install ${missing_packages[*]}"
            echo ""
            echo "Or using conda:"
            echo "conda install -y ${missing_packages[*]}"
            echo ""
        fi
    fi
}

# 创建环境配置文件
create_environment_config() {
    local config_file="hic_pipeline_env.conf"
    
    blue "\nCreating environment config file: $config_file"
    
    cat > "$config_file" << EOF
# Hi-C Pipeline environment config file
# Generated at: $(date)

# Software paths
PYTHON3_PATH=$(command -v python3 2>/dev/null || echo "/usr/bin/python3")
SAMTOOLS_PATH=$(command -v samtools 2>/dev/null || echo "/usr/bin/samtools")
BOWTIE2_PATH=$(command -v bowtie2 2>/dev/null || echo "/usr/bin/bowtie2")
FASTQC_PATH=$(command -v fastqc 2>/dev/null || echo "/usr/bin/fastqc")

# Optional tool paths
COOLER_PATH=$(command -v cooler 2>/dev/null || echo "")
MULTIQC_PATH=$(command -v multiqc 2>/dev/null || echo "")

# Python environment
PYTHON_VERSION=$(python3 --version 2>/dev/null || echo "unknown")
PYTHON_PACKAGES="$(pip3 list 2>/dev/null | grep -E '(pysam|bx-python|numpy|pandas|matplotlib)' | tr '\n' ';')"

# System info
SYSTEM_INFO="$(uname -a)"
CPU_INFO="$(grep 'model name' /proc/cpuinfo 2>/dev/null | head -n1 | cut -d':' -f2 | sed 's/^[[:space:]]*//')"
MEMORY_INFO="$(free -h 2>/dev/null | grep '^Mem:' | awk '{print \$2}' || echo 'unknown')"

# Recommended settings
RECOMMENDED_CPU=4
RECOMMENDED_MEMORY="8G"
RECOMMENDED_DISK_SPACE="50G"

# Check status
CHECK_STATUS="$(date): Environment check completed"
EOF
    
    green "Environment config file created: $config_file"
}

# 主检查函数
main_check() {
    echo "========================================"
    blue "Hi-C Pipeline environment check"
    blue "Based on nf-core/hic pipeline"
    echo "========================================"
    
    local all_good=true
    
    # 检查必需命令
    blue "\nChecking required command-line tools:"
    for cmd in "${REQUIRED_COMMANDS[@]}"; do
        if ! check_command "$cmd"; then
            all_good=false
        fi
