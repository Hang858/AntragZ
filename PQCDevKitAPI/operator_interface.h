/**
 * 更新时间： 2026-01-12 14:31
 * @file operator_interface.h
 * @brief 硬件算子抽象接口汇总，用于软件团队适配算法工作。
 * * 依据《硬件算子汇总》文档整理。
 */

#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

 /* --- 宏定义与常量 --- */

//待补充
#define OP_SUCCESS 0
#define OP_FAILURE -1

#define MAX_NTT_OMEGA_COUNT 8  // 最大支持8个旋转因子（4个正变换，4个逆变换）

//1:SHA256; 2:SM3  3:SHAKE256; 待定：4:SHAKE128; 5:SHA3-256; 6: SHA3-384; 7: SHA3-512
#define OP_ALG_SHA256   1
#define OP_ALG_SM3      2
#define OP_ALG_SHAKE128 4
#define OP_ALG_SHAKE256 3
#define OP_ALG_SHA3_256 5
#define OP_ALG_SHA3_384 6
#define OP_ALG_SHA3_512 7

#define OP_MODE_NORMAL 0
#define OP_MODE_LMS 1
#define OP_MODE_SPHINCS 2

#define ABORT_IF_FAIL(x) if ((x) != OP_SUCCESS) { \
        abort(); \
    } \

 /* --- 硬件算子接口声明 --- */
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
int OP_ntt256_init(const int32_t *omega, int32_t count, int32_t q, int32_t inv_q, int32_t inv_n, int32_t k);


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
int OP_ntt256(int32_t a_out[256], const int32_t a_in[256], int omega_index);


/**
 * @brief 算子1.2. 逆数论变换 (INTT)  
 * * 将多项式从频域转换回时域。
 * 
 * @param a_out  输出：恢复后的时域多项式 (256个 32-bit 元素)
 * @param a_in   输入：频域多项式系数 (256个 32-bit 元素)
 * @param omega_index 输入：选择使用的旋转因子索引 (count/2 到 count - 1)
 * @return       0: 成功; -1: 失败。
 */
int OP_intt256(int32_t a_out[256], const int32_t a_in[256], int omega_index);


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
int OP_cwm(int32_t *c_out, const int32_t *a_in, const int32_t *b_in, int32_t n, int32_t q, int32_t k);



/**
 * @brief 算子2. 向量乘法  (元素数量可变)
 * * Z = (X * Y) mod q
 * 
 * @param z_out     输出：一个 16-bit 元素。
 * @param x_in      输入：向量  (length个 16-bit 元素)。
 * @param y_in      输入：向量  (length个 16-bit 元素)。
 * @param length    输入：向量长度。（最大长度：Frodo: 1344, scloudplus: 1136）
 * @param q         输入：模数 (16-bit) 2的次幂。
 * @return          0: 成功; -1: 失败。
 * 待定：数据无符号
 */
int OP_vector_mul(uint16_t *z_out, const uint16_t *x_in, const uint16_t *y_in, uint16_t length, uint16_t q);


/**
 * @brief 算子3.1 单包哈希
 * 
 * @param alg        输入：决定使用的哈希算法类型。(1:SHA256; 2:SM3  3:SHAKE256; 待定：4:SHAKE128; 5:SHA3-256; 6: SHA3-384; 7: SHA3-512)
 * 
 * @param mode       输入：决定使用的哈希算法模式。(0:普通； 1:哈希链)
 * 
 * 注意：哈希链模式下只支持 SM3 和 SHAKE256
 * 
 * @param n          输入：输出哈希值的长度 (字节)。(SHAKE算法支持任意字节，算法待定，例如：SHAKE256最大输出1000字节，其他算法不大于算法固定输出长度)
 * @param input      输入：待哈希的数据缓冲区指针。
 * @param input_len  输入：待哈希数据的长度 (字节)。
 * @param link_count    输入：哈希链模式下哈希次数。(仅哈希链使用，其他默认0)
 * @param output     输出：存储哈希结果的缓冲区指针。
 * @return           0: 成功; -1: 失败。
 */
int OP_hash(uint8_t alg, uint8_t mode, int n, const void *input, int input_len, uint8_t link_count, void *output);

/**
 * @brief 算子3.2 多包哈希初始化（软件实现：清0）
 * 
 * @param alg        输入：决定使用的哈希算法类型。(1:SHA256; 2:SM3  3:SHAKE256; 待定：4:SHAKE128; 5:SHA3-256; 6: SHA3-384; 7: SHA3-256)
 * @param s      输入/输出：哈希的状态。
 * @param s_len      输入：哈希的状态的长度 (字节)。
 * @return           0: 成功; -1: 失败。
 */
int OP_hash_init(uint8_t alg, void *s, int s_len);
/**
 * @brief 算子3.3 多包哈希海绵吸收函数（可以多次调用）
 * 
 * @param alg        输入：决定使用的哈希算法类型。(1:SHA256; 2:SM3  3:SHAKE256; 待定：4:SHAKE128; 5:SHA3-256; 6: SHA3-384; 7: SHA3-256)
 * @param s      输入/输出：哈希的状态。
 * @param s_len      输入：哈希的状态的长度 (字节)。
 * @param input      输入：待哈希的数据缓冲区指针。
 * @param input_len  输入：待哈希数据的长度 (字节)。
 * @return           0: 成功; -1: 失败。
  */
int OP_hash_absorb(uint8_t alg, void *s, int s_len, const void *input, int input_len);
/**
 * @brief 算子3.4 多包哈希输出函数（SHAKE可以多次调用，其他算法只能调用一次）
 * 
 * @param alg        输入：决定使用的哈希算法类型。(1:SHA256; 2:SM3  3:SHAKE256; 待定：4:SHAKE128; 5:SHA3-256; 6: SHA3-384; 7: SHA3-256)
 * @param s      输入/输出：哈希的状态。
 * @param s_len      输入：哈希的状态的长度 (字节)。
 * @param output     输出：存储哈希结果的缓冲区指针。
 * @param output_len      输入：存储哈希结果的缓冲区的长度 (字节)。
 * @return           0: 成功; -1: 失败。
  */
int OP_hash_squeeze(uint8_t alg, void *s, int s_len, void *output, int output_len);



 /**
 * @brief 算子4. 硬件随机数生成  (拒绝均匀采样)
 * 内部实现：
 * 调用硬件TRNG模块，生成真随机数作为种子；
 * 基于该种子，生成伪随机数。
 * 
 * @param buffer  输出：存储随机数的缓冲区。
 * @param size    输入：需要获取的随机数长度。
 * @return        0: 成功; -1: 硬件故障/熵源不足。
 * 待定：接口是否直接调用硬件TRNG模块？
 * 待定：试点算法未定
 */
int OP_trng(void *buffer, int size);


/**
 * @brief 算子5. 有限域运算
 *  0: 模加        C=(a+b) mod n
 *  1: 模减        C=(a-b) mod n
 *  2: 模乘        C=(a*b) mod n
 *  3: 模幂        C=a^b   mod n
 *  4: 模逆        C=a^(-1)mod n
 * @param c_out    输出：结果 C (大整数，数组表示，小端字节序)
 * @param opr 输入：运算类型 (0: 模加；1: 模减；2: 模乘；3: 模幂；4：模逆)
 * @param a_in     输入：底数 a (大整数，数组表示，小端字节序)
 * @param b_in     输入：指数 b (大整数，数组表示，小端字节序)
 * @param n_in     输入：模数 n (大整数，数组表示，小端字节序，软件实现仅需支持5*2^248-1、65*2^376-1、27*2^500-1)
 * @param len      输入：大整数的长度，每个数组有多少个 64-bit 元素
 * @return         0: 成功; -1: 失败。
 */
int OP_finite_field(uint64_t *c_out, uint8_t opr, const uint64_t *a_in, const uint64_t *b_in, const uint64_t *n_in, int len);