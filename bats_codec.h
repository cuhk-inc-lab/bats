#ifndef BATS_CODEC_H
#define BATS_CODEC_H

/*
 * 给别的程序调用的接口。不包含丢包、打印和四节点仿真。
 *
 * 一个 batch 发出 BATS_M 个编码包。每个包是 BATS_M 字节系数加上 T 字节载荷。
 * 系数按包主序排：第 j 个包从 coeff[j * BATS_M] 开始。
 * 载荷同样按包主序：第 j 个包从 payload[j * T] 开始。
 * G 不在这些字节里。源和目的用同一个 code_seed 加上 batch_id 各自复原。
 *
 * 中继一次只处理一个 batch。不同 batch 的包不要放进同一次 bats_recode。
 * 调用方自己保存 recode_state，同一跳上的后续 batch 继续用这个状态。
 *
 * 编译进调用方：
 *   gcc -std=c11 -c field.c source.c dest.c codec.c
 *   ar rcs libbats.a field.o source.o dest.o codec.o
 * 调用方只包含本头文件。使用前调用一次 bats_init。
 */

#include <stddef.h>
#include <stdint.h>

#define BATS_M 8

typedef struct {
    uint32_t batch_id;
    uint8_t coeff[BATS_M];
    const uint8_t *payload; /* T 字节，只在 bats_decode 返回前有效 */
} BatsSymbol;

void bats_init(void);

/* 写出 BATS_M 个包。成功返回 0，参数不合法返回 -1。 */
int bats_encode(const uint8_t *src, int K, int T, uint64_t code_seed, uint32_t batch_id,
                uint8_t *coeff, uint8_t *payload);

/*
 * 用已经收到的 n_in 个包混合出新包。n_in 为 0 时不写出包，*n_out 为 0。
 * 否则写出 BATS_M 个包，*n_out 为 BATS_M。
 * coeff 与 payload 可以和输入重叠。成功返回 0，参数不合法返回 -1。
 */
int bats_recode(const uint8_t *coeff_in, const uint8_t *payload_in, int n_in, int T,
                uint64_t *recode_state, uint8_t *coeff_out, uint8_t *payload_out,
                int *n_out);

/*
 * symbols 可以来自多个 batch，按 batch_id 分组后再译码。
 * 解出全部 K 个源包返回 1，解不出返回 0，参数不合法返回 -1。
 * n_bp 和 n_inact 可以为 NULL。
 */
int bats_decode(const BatsSymbol *symbols, int n_symbols, int K, int T, uint64_t code_seed,
                uint8_t *dst, int *n_bp, int *n_inact);

#endif
