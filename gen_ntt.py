# gen_ntt_correct.py
def power(a, b, m):
    res = 1
    a %= m
    while b > 0:
        if b % 2 == 1:
            res = (res * a) % m
        a = (a * a) % m
        b //= 2
    return res

def modInverse(n, m):
    return power(n, m - 2, m)

def bit_reverse(n, bits):
    result = 0
    for i in range(bits):
        result = (result << 1) | (n & 1)
        n >>= 1
    return result

def main():
    Q = 268435457      # 0x10000001
    
    # 1. 计算关键常数
    # R = 2^32 mod Q
    # 2^32 = 16 * Q - 16  => R = -16 = 268435441
    R = 4294967296 % Q 
    
    # Q_inv mod 2^32
    # Q * Q_inv = 1 mod 2^32
    # 0x10000001 * 0xF0000001 = ...0001
    # Q_inv = 4026531841
    # 验证: (Q * 4026531841) % 4294967296 should be 1
    Q_inv = pow(Q, -1, 4294967296)

    print(f"// Correct Constants for Q = {Q}")
    print(f"#define BIG_Q         {Q}")
    print(f"#define BIG_INV_Q     {Q_inv}U  // 0x{Q_inv:08X}")
    print(f"#define R_MOD_Q       {R}       // -16")
    print(f"#define R_SQUARED_MOD_Q {(R*R)%Q}      // 256")
    
    # BIG_INV_2_MONT = (1/2) * R mod Q
    inv_2 = modInverse(2, Q)
    inv_2_mont = (inv_2 * R) % Q
    print(f"#define BIG_INV_2_MONT {inv_2_mont}  // -8 mod Q")
    print("-" * 40)

    # 2. 生成表
    N_EXT = 512
    # 原根 g=3
    g = 3
    w_512 = power(g, (Q - 1) // N_EXT, Q)
    inv_w_512 = modInverse(w_512, Q)
    w_256 = power(w_512, 2, Q)
    inv_w_256 = modInverse(w_256, Q)

    # BIG_OMEGA_TABLE_256 (BitReversed, Montgomery Form)
    fwd_table = [(power(w_256, bit_reverse(i, 8), Q) * R) % Q for i in range(256)]
    inv_table = [(power(inv_w_256, bit_reverse(i, 8), Q) * R) % Q for i in range(256)]

    print("static const int32_t BIG_OMEGA_TABLE_256[] = {")
    print("    // --- Forward (Mont) ---")
    for i in range(0, 256, 8):
        print("    " + ", ".join(f"{x}" for x in fwd_table[i:i+8]) + ",")
    print("    // --- Inverse (Mont) ---")
    for i in range(0, 256, 8):
        print("    " + ", ".join(f"{x}" for x in inv_table[i:i+8]) + ",")
    print("};")

    # PHI_512_FWD (Natural, Montgomery Form)
    phi_fwd = [(power(w_512, i, Q) * R) % Q for i in range(256)]
    print("\nstatic const int32_t PHI_512_FWD[256] = {")
    for i in range(0, 256, 8):
        print("    " + ", ".join(f"{x}" for x in phi_fwd[i:i+8]) + ",")
    print("};")

    # PHI_512_INV (Natural, Montgomery Form)
    phi_inv = [(power(inv_w_512, i, Q) * R) % Q for i in range(256)]
    print("\nstatic const int32_t PHI_512_INV[256] = {")
    for i in range(0, 256, 8):
        print("    " + ", ".join(f"{x}" for x in phi_inv[i:i+8]) + ",")
    print("};")

if __name__ == "__main__":
    main()