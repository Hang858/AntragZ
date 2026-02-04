import math
import sys

# ================= 配置区域 =================

# 1. 基础参数 (来自设计文档 Table 1)
# 这是基础采样器的标准差 r0。
# 所有的更大宽度的采样 (s0, Lr0) 都将由这个基础分布合成。
SIGMA = 1.27

# 2. 采样分解基数 (Beta)
# Zalcon 论文参数 Beta=64 (2^6)
BETA = 64

# 3. 分解层数 (Layers)
# Zalcon 论文参数 l=5
# 这意味着随机化舍入的精度为 Beta^l = 64^5 = 2^30 (30-bit)
LAYERS = 5

# 4. 概率精度 (64位)
# 使用 uint64_t 存储 CDF，满足工程安全性需求
BIT_PRECISION = 64
SCALE_FACTOR = 2 ** BIT_PRECISION
MAX_VAL = SCALE_FACTOR - 1

# 5. 截断尾部 (Tail Cut)
# Sigma=1.27 时，20*Sigma 已经非常安全
TAIL_CUT = 20

# 6. 输出文件名
OUTPUT_FILE = "cdt_zalcon_data.h"

# ================= 计算逻辑 =================

def gaussian_pdf(x, sigma, mu):
    """
    计算离散高斯分布的概率质量
    """
    return math.exp(-((x - mu) ** 2) / (2 * (sigma * sigma)))

def generate_single_table(mu):
    """
    为特定的中心 mu 生成单列 CDF 表
    """
    # 统一使用 [-TAIL_CUT, TAIL_CUT] 的相对索引
    min_val = -TAIL_CUT
    max_val = TAIL_CUT
    
    integers = range(min_val, max_val + 1)
    probs = []
    total_mass = 0.0
    
    # 1. 计算 PDF
    for x in integers:
        p = gaussian_pdf(x, SIGMA, mu)
        probs.append(p)
        total_mass += p
        
    # 2. 计算 CDF 并缩放
    cdf_ints = []
    current_cdf = 0.0
    
    for p in probs:
        current_cdf += p / total_mass
        scaled = int(current_cdf * SCALE_FACTOR)
        if scaled > MAX_VAL:
            scaled = MAX_VAL
        cdf_ints.append(scaled)
    
    # 强制修正最后一个值为最大值 (CDF=1.0)
    cdf_ints[-1] = MAX_VAL
    
    return cdf_ints

def generate_all_cdts():
    """
    生成所有需要的 CDT 表
    根据对称性，只需要生成 mu = 0/BETA 到 (BETA/2)/BETA 的表
    即 [0, 1/64, ..., 32/64]
    """
    print(f"[Info] Configuration: Sigma={SIGMA}, Beta={BETA}, Layers={LAYERS}")
    
    tables = []
    # 利用对称性 D_(1-c) = 1 - D_c，只存储 [0, 0.5]
    num_tables = (BETA // 2) + 1
    
    print(f"[Info] Generating {num_tables} tables for centers 0.0 to 0.5...")
    
    for k in range(num_tables):
        mu = k / float(BETA)
        table = generate_single_table(mu)
        tables.append(table)
        
    return tables

def write_c_header(tables):
    """写入 C 头文件"""
    num_tables = len(tables)
    table_len = len(tables[0])
    min_val = -TAIL_CUT
    max_val = TAIL_CUT
    
    content = f"""/**
 * @file {OUTPUT_FILE}
 * @brief Auto-generated CDT tables for Zalcon/ANTRAG-Z
 * * Parameters:
 * - Sigma (r0): {SIGMA}
 * - Beta:       {BETA}
 * - Layers (l): {LAYERS} (Precision = Beta^l = 2^30)
 * - Precision:  {BIT_PRECISION}-bit (uint64_t)
 * - Range:      [{min_val}, {max_val}]
 * * Note: Stores CDFs for centers [0, 1/{BETA}, ..., 0.5].
 * Use symmetry D_(1-c) for c > 0.5.
 */

#ifndef CDT_ZALCON_DATA_H
#define CDT_ZALCON_DATA_H

#include <stdint.h>

#define CDT_SIGMA     {SIGMA}
#define CDT_BETA      {BETA}
#define CDT_LAYERS    {LAYERS}  // Zalcon parameter l
#define CDT_COLS      {num_tables}   // 0 to Beta/2
#define CDT_ROWS      {table_len}    // Size of CDF
#define CDT_OFFSET    {abs(min_val)} // Index offset
#define CDT_MIN_VAL   ({min_val})

/**
 * CDT Lookup Table
 * Data type: uint64_t (Scaled by 2^64)
 */
static const uint64_t CDT_TABLE[CDT_COLS][CDT_ROWS] = {{
"""
    
    for k, table in enumerate(tables):
        mu = k / float(BETA)
        content += f"    // Center mu = {k}/{BETA} ({mu:.4f})\n"
        content += "    {\n"
        
        chunk_lines = []
        for i in range(0, len(table), 4):
            chunk = table[i:i+4]
            hex_strs = [f"0x{x:016X}ULL" for x in chunk] 
            chunk_lines.append("        " + ", ".join(hex_strs))
            
        content += ",\n".join(chunk_lines)
        
        if k < num_tables - 1:
            content += "\n    },\n"
        else:
            content += "\n    }\n"
            
    content += "};\n\n#endif // CDT_ZALCON_DATA_H\n"
    
    with open(OUTPUT_FILE, "w") as f:
        f.write(content)
    
    print(f"[Success] Generated {OUTPUT_FILE}")

# ================= 主程序 =================

if __name__ == "__main__":
    all_tables = generate_all_cdts()
    write_c_header(all_tables)