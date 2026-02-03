#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "operator_interface.h"

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
int OP_trng(void *buffer, int size)
{
    srand((unsigned int)time(NULL)); // 使用当前时间作为随机数种子
    unsigned char *buf = (unsigned char *)buffer;
    for (int i = 0; i < size; i++) {
        buf[i] = rand() % 256; // 生成0-255之间的随机字节
    }
    return OP_SUCCESS;
}
