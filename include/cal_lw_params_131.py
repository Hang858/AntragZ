import math

# ================= 配置参数 =================
r0 = 1.27
target_s = 131.0

# ================= 辅助函数 =================

def get_next_r(r, z):
    """计算下一层的标准差: r_new = r * sqrt(z^2 + (z-1)^2)"""
    return r * math.sqrt(z**2 + (z-1)**2)

def find_best_sequence(target, max_depth=3):
    """搜索最优的 z 序列"""
    best_err = float('inf')
    best_seq = []
    best_r = 0
    
    # 深度优先搜索所有可能的 z 组合
    # 限制 z 的范围以加速搜索，并确保增长合理
    
    # Depth 1 (一般达不到 131，但作为基准)
    for z1 in range(2, 200):
        r1 = get_next_r(r0, z1)
        err = abs(r1 - target)
        if err < best_err:
            best_err = err
            best_seq = [z1]
            best_r = r1
            
    # Depth 2
    for z1 in range(2, 50):
        r1 = get_next_r(r0, z1)
        if r1 > target * 1.5: continue
        
        # 估算 z2
        z2_approx = (target / r1) / math.sqrt(2)
        z2_center = max(2, round(z2_approx))
        
        for z2 in range(z2_center - 2, z2_center + 3):
            if z2 < 2: continue
            r2 = get_next_r(r1, z2)
            err = abs(r2 - target)
            if err < best_err:
                best_err = err
                best_seq = [z1, z2]
                best_r = r2

    # Depth 3 (你的目标可能需要3层)
    for z1 in range(2, 20):
        r1 = get_next_r(r0, z1)
        for z2 in range(2, 20):
            r2 = get_next_r(r1, z2)
            if r2 > target * 1.5: continue
            
            z3_approx = (target / r2) / math.sqrt(2)
            z3_center = max(2, round(z3_approx))
            
            for z3 in range(z3_center - 2, z3_center + 3):
                if z3 < 2: continue
                r3 = get_next_r(r2, z3)
                err = abs(r3 - target)
                if err < best_err:
                    best_err = err
                    best_seq = [z1, z2, z3]
                    best_r = r3

    return best_seq, best_r, best_err

def generate_weights_and_coeffs(z_seq):
    """根据 z 序列生成权重和系数"""
    weights = [1]
    coeffs = [1]
    
    for z in z_seq:
        new_weights = []
        new_coeffs = []
        for w, a in zip(weights, coeffs):
            # 权重分裂: w -> w*z, w*(z-1)
            new_weights.append(w * z)
            new_weights.append(w * (z - 1))
            
            # 系数分裂 (Bézout系数): 
            # 使得 a'*(w*z) + b'*(w*(z-1)) = a*w
            # 解为: a' = a, b' = -a (这是一个特解，且系数绝对值最小为1)
            # 验证: a*z + (-a)*(z-1) = az - az + a = a
            new_coeffs.append(a)
            new_coeffs.append(-a)
            
        weights = new_weights
        coeffs = new_coeffs
    return weights, coeffs

# ================= 主流程 =================

print(f"=== SampleSW 参数计算 (Target s={target_s}) ===")
seq, r_final, err = find_best_sequence(target_s, 3)

print(f"最优 z 序列: {seq}")
print(f"实际标准差: {r_final:.4f}")
print(f"误差: {err:.4e} ({err/target_s:.2%})")

weights, coeffs = generate_weights_and_coeffs(seq)

print("\n=== C Language Output ===")
print(f"// Generated for s={target_s}, r0={r0}")
print(f"// Sequence z: {seq}")
print(f"#define SW_LAYERS {len(seq)}")
print(f"#define SW_NUM_SAMPLES {1 << len(seq)}")

print("static const int64_t SW_WEIGHTS[SW_NUM_SAMPLES] = {")
print("    " + ", ".join(f"{w}LL" for w in weights))
print("};")

print("static const int64_t SW_COEFFS[SW_NUM_SAMPLES] = {")
print("    " + ", ".join(f"{c}LL" for c in coeffs))
print("};")