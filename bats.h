#ifndef BATS_H
#define BATS_H

/*
 * BATS 线网仿真（只在本进程内存里传）
 *
 * 四节点：源 → 中继1 → 中继2 → 目的。
 * 无预编码、无 socket、无端口。丢包是入队前按概率丢掉。
 *
 * 外码：每个 batch 从写死的 Ψ 抽度数 d，均匀选 d 个源包，
 *       乘 d×M 的随机矩阵 G。G 不进包。
 *       目的端用同一个种子和 batch_id 自己算出 d、下标和 G。
 * 第 j 个发出包的系数向量是 I_M 的第 j 列。payload 已是 (BG) 的第 j 列，
 * 不再乘一次单位阵。
 * 中继：只混合同一个 batch_id 里已经收到的包（没收齐也 recode），
 *       payload 和系数乘同一个随机矩阵。不同 batch 不混。
 * 目的：不 recode。先 BP（rank 等于当前度数就解），解不动再 inactivation。
 *
 * 按节点分文件：
 *   source.c  源。抽度数、生成 G、发出 M 个包
 *   relay.c   中继。中继1 和中继2 跑同一段 recode
 *   dest.c    目的。收包、BP、inactivation
 *   sim.c     把上面三个节点串成 源 → 中继1 → 中继2 → 目的
 *   link.c    节点之间的内存链路
 *   field.c   GF(256) 和随机数，三个节点共用
 *   main.c    自检、例子、开销测量、文件还原
 *
 * 用法：
 *   make && ./bats_sim [输入文件] [输出文件]
 *   build.bat              （Windows，用本目录里的 tcc）
 * 还原出的文件与原文一致时返回 0，否则返回 1。
 * 另外对 K=128 从少到多搜索刚好能还原的 batch 数，并打印每跳秩分布和相对 K/h_bar 的开销。
 */

#define _CRT_SECURE_NO_WARNINGS

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define M 8
#define LOSS_PERCENT 5
#define FILE_T 32

void *xmalloc(size_t n);
void *xrealloc(void *p, size_t n);

void gf_init(void);
uint8_t gf_mul(uint8_t a, uint8_t b);
uint8_t gf_inv(uint8_t a);
int gf_self_test(void);

/* A 是 nU×nR。成功时 Inv * S = I，S 的第 i 列是 A 的第 col_of_pivot[i] 列。返回秩。 */
int row_reduce(const uint8_t *A, int nU, int nR, int *col_of_pivot, uint8_t *Inv);

typedef struct {
    uint64_t s;
} Rng;

uint32_t rng_below(Rng *r, uint32_t n);
uint8_t rng_byte(Rng *r);
Rng rng_for_batch(uint64_t seed, uint32_t batch_id);

extern const int PSI_DEG[M];
extern const int PSI_W[M];
extern const int PSI_SUM;

typedef struct {
    int d;
    int *sel;
    uint8_t *G; /* d 行 M 列，行主序 */
} Plan;

Plan make_plan(uint64_t seed, uint32_t batch_id, int K);
void plan_free(Plan *plan);
int batches_needed(int K);

typedef struct {
    uint32_t batch_id;
    uint8_t coeff[M];
    uint8_t *payload;
} Packet;

typedef struct {
    Packet *p;
    int n;
    int cap;
} Vec;

typedef struct {
    int sent;
    int dropped;
    int delivered;
} LinkStat;

void vec_push(Vec *v, Packet pkt);
void vec_free(Vec *v);
void print_coeff_payload(const char *tag, uint32_t batch_id, int j,
                         const uint8_t coeff[M], const uint8_t *payload,
                         int T, const char *note);
int link_send(Packet pkt, Vec *dst, Rng *loss, int loss_percent, LinkStat *st,
              int do_print, const char *link, int j, int T);

/* 收到的系数向量排成 r×M，秩就是这一跳该 batch 的 transfer 秩。整批丢失为 0。 */
int coeff_rank(const Packet *pkts, int r);

/* 源：这个 batch 编码出 M 个包，经链路交给下一跳。 */
void source_send(const uint8_t *src, int K, int T, uint64_t code_seed, int batch,
                 Vec *next, Rng *loss, int loss_percent, LinkStat *st,
                 int verbose, int trace, const char *link_name);

/*
 * 中继：只混合同一个 batch。空批不 recode。
 * mention_m 为真时，收到包数的打印里带上 M（中继1 的原文如此）。
 * 有包时会释放 in。
 */
void relay_forward(Vec *in, Vec *next, Rng *recode, Rng *loss, int loss_percent,
                   LinkStat *st, int batch, int T, int verbose, int trace,
                   const char *name, const char *link_name, int mention_m,
                   int *rank_out);

typedef struct {
    int *nb;
    uint8_t *gamma; /* n_nb * n_recv */
    uint8_t *Y;     /* n_recv * T */
    int n_nb;
    int cap_nb;
    int n_recv;
    int batch_id;
} BatchEq;

void batch_free(BatchEq *b);
void eqs_from_packets(BatchEq *eq, const Vec *got, uint64_t code_seed, int batch_id, int K, int T);

/* 目的：不 recode。把收到的包收成方程，并释放 got。r==0 时不释放（本来就是空的）。 */
void dest_collect(BatchEq *eq, Vec *got, uint64_t code_seed, int batch, int K, int T,
                  int verbose, int trace, int *rank_out);

int decode(BatchEq *batches, int n_batches, int K, int T, uint8_t *dst,
           int *n_bp_out, int *n_inact_out);
int inactivation_self_test(void);

typedef struct {
    int ok;
    int bp_solved;
    int inactivated;
    LinkStat link[3];
} SimResult;

SimResult simulate_ex(const uint8_t *src, int K, int T, uint64_t code_seed,
                      uint64_t loss_seed, uint64_t recode_seed, int loss_percent,
                      int n_batches, int verbose, int trace_batch, uint8_t *dst,
                      int *rank_r1, int *rank_r2, int *rank_dst, Vec **dest_keep);
SimResult simulate(const uint8_t *src, int K, int T, uint64_t code_seed,
                   uint64_t loss_seed, uint64_t recode_seed, int loss_percent,
                   int n_batches, int verbose, int trace_batch, uint8_t *dst);

#endif
