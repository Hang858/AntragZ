#include "operator_interface.h"


// 全局变量和状态可以在这里定义
int32_t g_ntt_omega[MAX_NTT_OMEGA_COUNT * 256]={0}; // 存储逆旋转因子
int32_t g_ntt_count=0;      // 存储逆旋转因子数组长度
int32_t g_ntt_q=0;          // 存储模数
int32_t g_ntt_inv_q=0;
int32_t g_ntt_inv_n=0;
int32_t g_ntt_k=0;          // 存储剩余层数


 /* --- 硬件算子接口实现 --- */
/**
 * @brief 算子1.0. 初始化
 * @param omega  输入：旋转因子数组 (count个旋转因子，每个旋转因子包含256个 32-bit 元素，前面是所有的正变换旋转因子，后面是所有的逆变换旋转因子)
 * @param count 输入：旋转因子的个数（2的倍数，正变换和逆变换使用的旋转因子不同，最多支持8个旋转因子）
 * @param q      输入：模数(参考8380417)
 * @param inv_q  输入：q^(-1)模2^32
 * @param inv_n  输入：256^(-1)模q
 * @param k      输入：剩余层数（0,, 1, 2, 3）
 * @return       0: 成功; -1: 失败。
 */
int OP_ntt256_init(const int32_t *omega, int32_t count, int32_t q, int32_t inv_q, int32_t inv_n, int32_t k)
{
    int i;  
    if (omega==NULL||count<=0||count>MAX_NTT_OMEGA_COUNT||k<0||k>3) return OP_FAILURE;
    for(i=0;i<256*count;i++) g_ntt_omega[i]=omega[i];
    g_ntt_count = count;
    g_ntt_q = q;
    g_ntt_inv_n = inv_n;
    g_ntt_inv_q = inv_q;
    g_ntt_k = k;
    return OP_SUCCESS;
}

//补充函数：Montgomery映射函数，输入64bit乘积a，输出32bit整数r，满足r = a * Beta^{-1} mod q，其中Beta=2^{32}
static int32_t MontgomeryReduce(int64_t a, int32_t q, int32_t inv_q)  
{
    int64_t t;
    int64_t m,r;
    int64_t beta = 1ULL << 32; 
    m = a * inv_q;
    m=m&(beta -1);  // m mod β  
    t = m * q;
    r = a - t;
    r >>= 32; // r in thr range of (-Q, Q)
    r = r + ((r >> 31) & q); // Ensure r is in the range of (0, Q)
    return (int32_t)r;
}
 /**
 * @brief 算子1.1. 正向数论变换 (NTT)  
 * * 将多项式从时域转换到频域。
 * * 待定：分解到倒数第三层，或者增加参数动态确定停止层数。
 * 
 * @param a_out  输出：变换后的系数数组 (256个 32-bit 元素)
 * @param a_in   输入：待变换的多项式系数数组 (256个 32-bit 元素)
 * @param omega_index 输入：选择使用的旋转因子索引 (0 到 count/2 - 1)
 * @return       0: 成功; -1: 失败。
 */
int OP_ntt256(int32_t a_out[256], const int32_t a_in[256], int omega_index)
{
    if(a_in == NULL)
    {
        return -1; //错误处理：输入指针为空
    }

    int t,m,i,j,s,e;
    int32_t U=0,V=0,S;
    t=256; 
    int end_level = 256 >> g_ntt_k;
    int start_index = omega_index<<(8-g_ntt_k);
    for(i=0;i<256;i++)
        a_out[i]=a_in[i];
    for(m=1;m<end_level;m<<=1) 
    {
        t=(t>>1);
        for(i=0;i<m;i++)  
        {
            s=(i*t)<<1;
            e=s+t;
            S=g_ntt_omega[start_index +m+i]; //单位根通过初始化操作预先存储在omega数组中
            for(j=s;j<e;j++)  
            {
                U=a_out[j];
                V=a_out[j+t];

                V=MontgomeryReduce((int64_t)V*(int64_t)S, g_ntt_q, g_ntt_inv_q); 

                //Map each coefficient to its standard representation.
                a_out[j]=U+V-g_ntt_q;
                a_out[j] = a_out[j] + ((a_out[j]>>31)&g_ntt_q);

                a_out[j+t]=U-V;
                a_out[j+t] = a_out[j+t] + ((a_out[j+t]>>31)&g_ntt_q);
            }
        }
    }

    return 0;
}


/**
 * @brief 算子1.2. 逆数论变换 (INTT)  
 * * 将多项式从频域转换回时域。
 * 
 * @param a_out  输出：恢复后的时域多项式 (256个 32-bit 元素)
 * @param a_in   输入：频域多项式系数 (256个 32-bit 元素)
 * @param omega_index 输入：选择使用的旋转因子索引 (count/2 到 count - 1)
 * @return       0: 成功; -1: 失败。
 */
int OP_intt256(int32_t a_out[256], const int32_t a_in[256], int omega_index)
{
    if(a_in == NULL)
    {
        return -1; //错误处理：输入指针为空
    }

    int t,m,i,j,s,e,h;
    int32_t U=0,V=0,S;

    for(i=0;i<256;i++)
        a_out[i]=a_in[i];
    
    t = 1 << g_ntt_k;
    int begin_level = 256 >> g_ntt_k;
    int start_index = omega_index*(256>>g_ntt_k);

    for (m = begin_level; m > 1; m >>= 1)
    {
        s=0;
        h=m>>1;
        for(i=0;i<h;i++) 
        {
            e=s+t;
            S=g_ntt_omega[start_index + h+i]; //单位根通过初始化操作预先存储在omega数组中
            for(j=s;j<e;j++)  
            {
                U=a_out[j];
                V=a_out[j+t];

                 // Update the first element of the butterfly
                a_out[j] = (U + V - g_ntt_q);
                a_out[j] = a_out[j] + ((a_out[j] >> 31) & g_ntt_q); // Ensure the result is in the range (0, Q)

                // Update the second element of the butterfly using Montgomery reduction
                a_out[j + t] = MontgomeryReduce((int64_t)(U - V) * (int64_t)S, g_ntt_q, g_ntt_inv_q);
            }
            s=s+2*t;
        }
        t*=2;
    }
    
    for(i=0;i<256;i++) //归一化处理
    {
        a_out[i] = MontgomeryReduce((int64_t)a_out[i]* (int64_t)g_ntt_inv_n, g_ntt_q, g_ntt_inv_q); // Final Montgomery reduction to ensure standard representation
    }


    return 0;
}

/**
 * @brief 算子1.3. 多项式乘 (CWM - Coefficient-Wise Multiplication)
 * * 在频域对两个多项式进行对应系数模乘。
 * 
 * @param c_out  输出：点乘结果多项式  (n个 32-bit 元素)
 * @param a_in   输入：多项式 A 的频域表示 (n个 32-bit 元素)
 * @param b_in   输入：多项式 B 的频域表示 (n个 32-bit 元素)
 * @param n      输入：维数(最大 1024)
 * @param q      输入：模数(最长 31bit)
 * @param k      输入：剩余层数（硬件支持0, 1; 软件实现2, 3）
 * @return       0: 成功; -1: 失败。
 */
int OP_cwm(int32_t *c_out, const int32_t *a_in, const int32_t *b_in, int32_t n, int32_t q, int32_t k)
{
    int i, j;
    if (k == 0) {
        for (i = 0; i < n; i++)
        {
            c_out[i] = MontgomeryReduce((int64_t)a_in[i] * (int64_t)b_in[i], q, g_ntt_inv_q);
        }
    }
    else if (k == 1) {
        for (i = 0, j = 0; i < n - 3; i += 4, j++)
        {
            const int64_t r = g_ntt_omega[(n >> 2) + j];

            /* ---------- 第一对 ---------- */
            int64_t a0 = a_in[i];
            int64_t a1 = a_in[i + 1];
            int64_t s0 = b_in[i];
            int64_t s1 = b_in[i + 1];

            int64_t t00 = a0 * s0;                     // < Q^2
            int64_t t11 = a1 * s1;                     // < Q^2
            int64_t t11r = MontgomeryReduce(t11, q, g_ntt_inv_q) * r; // < Q^2

            c_out[i] = MontgomeryReduce(t00 + t11r, q, g_ntt_inv_q);

            int64_t t01 = a0 * s1;
            int64_t t10 = a1 * s0;

            c_out[i + 1] = MontgomeryReduce(t01 + t10, q, g_ntt_inv_q);

            /* ---------- 第二对 ---------- */
            int64_t a2 = a_in[i + 2];
            int64_t a3 = a_in[i + 3];
            int64_t s2 = b_in[i + 2];
            int64_t s3 = b_in[i + 3];

            int64_t u22 = a2 * s2;
            int64_t u33 = a3 * s3;
            int64_t u33r = MontgomeryReduce(u33, q, g_ntt_inv_q) * r;

            c_out[i + 2] = MontgomeryReduce(u22 - u33r, q, g_ntt_inv_q);

            int64_t u23 = a2 * s3;
            int64_t u32 = a3 * s2;

            c_out[i + 3] = MontgomeryReduce(u23 + u32, q, g_ntt_inv_q);
        }
    }
    // else if(k==2)
    //     basemul4t4(c_out, a_in, b_in, omega);
    // else if(k==3)
    //     basemul8t8(c_out, a_in, b_in, omega);
    else
        return -1; //错误处理：不支持的k值
    return 0;
}
