#ifndef BATS_INTERNAL_H
#define BATS_INTERNAL_H

/* 编码、recode、译码共用的内部声明。不含链路、丢包和四节点仿真。 */

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "bats_codec.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M BATS_M

void *xmalloc(size_t n);
void *xrealloc(void *p, size_t n);

void gf_init(void);
uint8_t gf_mul(uint8_t a, uint8_t b);
uint8_t gf_inv(uint8_t a);
int gf_self_test(void);

/* A 是 nU×nR。成功时 Inv * S = I，S 的第 i 列是 A 的第 col_of_pivot[i] 列。返回秩。 */
int row_reduce(const uint8_t *A, int nU, int nR, int *col_of_pivot, uint8_t *Inv);

/* 行主序的 rows×M 系数矩阵。rows<=0 时返回 0。 */
int coeff_rank_matrix(const uint8_t *coeff, int rows);

typedef struct {
    uint64_t s;
} Rng;

uint32_t rng_below(Rng *r, uint32_t n);
uint8_t rng_byte(Rng *r);
Rng rng_for_batch(uint64_t seed, uint32_t batch_id);

typedef struct {
    int d;
    int *sel;
    uint8_t *G; /* d 行 M 列，行主序 */
} Plan;

Plan make_plan(uint64_t seed, uint32_t batch_id, int K);
void plan_free(Plan *plan);
void encode_from_plan(const uint8_t *src, int T, const Plan *plan, uint8_t *coeff, uint8_t *payload);

typedef struct {
    int *nb;
    uint8_t *gamma; /* n_nb * n_recv */
    uint8_t *Y;     /* n_recv * T */
    int n_nb;
    int cap_nb;
    int n_recv;
    uint32_t batch_id;
} BatchEq;

void batch_free(BatchEq *b);

/* coeff、payload 按包主序，各 n_recv 个包。n_recv 为 0 时只记下 batch_id。 */
void batch_from_received(BatchEq *eq, const uint8_t *coeff, const uint8_t *payload, int n_recv,
                         uint64_t code_seed, uint32_t batch_id, int K, int T);

int decode(BatchEq *batches, int n_batches, int K, int T, uint8_t *dst,
           int *n_bp_out, int *n_inact_out);
int inactivation_self_test(void);

int psi_weight_sum(void);
long long psi_degree_moment(void);

#endif
