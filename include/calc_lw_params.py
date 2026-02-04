import math

# ================= 配置参数 =================
# 你的设计参数
r0 = 1.27
L_target_log2 = 35
L = 2 ** L_target_log2
Target_Width = r0 * L

# Zalcon 建议的前几层参数 (经验证这也适用于你的 r0)
# z 的选择原则是让 r 呈几何级数增长，避免某一层跳跃过大
# z0=3, z1=14, z2=304 是经过优化的序列
fixed_zs = [3, 14, 304] 

# ================= 计算逻辑 =================

def calc_next_r(r_prev, z):
    """计算下一层的标准差: r_new = r_prev * sqrt(z^2 + (z-1)^2)"""
    factor = math.sqrt(z**2 + (z - 1)**2)
    return r_prev * factor

def find_last_z(current_r, target_r):
    """
    已知当前 r 和目标 r，反推最后一层的 z
    r_target = current_r * sqrt(z^2 + (z-1)^2)
    近似: sqrt(2 * z^2) = ratio
    z * sqrt(2) = ratio
    z = ratio / sqrt(2)
    """
    ratio = target_r / current_r
    z_approx = ratio / math.sqrt(2)
    return round(z_approx)

def generate_weights(z_list):
    """
    根据 z 序列生成权重向量 v (长度 2^l)
    使用张量积逻辑: v_new = v_old (x) [z, z-1]
    """
    weights = [1]
    for z in z_list:
        new_weights = []
        for w in weights:
            # 每一层将之前的每个权重分裂为两个: w*z 和 w*(z-1)
            # 顺序: 偶数索引乘 z, 奇数索引乘 z-1 (符合二进制位为0和1的逻辑)
            new_weights.append(w * z)
            new_weights.append(w * (z - 1))
        weights = new_weights
    return weights

# ---------------- 主流程 ----------------

print(f"=== SampleLW 参数计算工具 ===")
print(f"基础标准差 r0: {r0}")
print(f"目标倍率 L: 2^{L_target_log2}")
print(f"目标宽度 r_target: {Target_Width:.4e}")
print("-" * 40)

# 1. 计算前 3 层
current_r = r0
all_zs = []

for i, z in enumerate(fixed_zs):
    all_zs.append(z)
    next_r = calc_next_r(current_r, z)
    print(f"Level {i+1} (z={z}): r = {next_r:.4f} (增长倍率: {next_r/current_r:.2f}x)")
    current_r = next_r

# 2. 计算最后一层 z3
z3 = find_last_z(current_r, Target_Width)
all_zs.append(int(z3))

# 3. 验证最终结果
final_r = calc_next_r(current_r, z3)
error_rate = (final_r - Target_Width) / Target_Width

print(f"Level 4 (z={z3}): r = {final_r:.4e}")
print("-" * 40)
print(f"最终 z 序列: {all_zs}")
print(f"目标宽度: {Target_Width:.6e}")
print(f"实际宽度: {final_r:.6e}")
print(f"相对误差: {error_rate:.8%}")

if abs(error_rate) > 0.01:
    print("[WARNING] 误差较大，请检查参数！")
else:
    print("[SUCCESS] 参数精度满足要求。")

# 4. 生成 C 语言数组
weights = generate_weights(all_zs)

print("\n=== C Language Code Output ===")
print(f"// Generated weights for SampleLW (L=2^{L_target_log2}, r0={r0})")
print(f"// Sequence z: {all_zs}")
print(f"#define LW_LAYERS 4")
print(f"#define LW_NUM_SAMPLES 16 // 2^4")
print("static const int64_t LW_WEIGHTS[LW_NUM_SAMPLES] = {")

# 格式化输出，每行4个
lines = []
for i in range(0, len(weights), 4):
    chunk = weights[i:i+4]
    line_str = ", ".join(f"{w}LL" for w in chunk)
    lines.append(f"    {line_str}")

print(",\n".join(lines))
print("};")