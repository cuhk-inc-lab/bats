#ifndef BATS_H
#define BATS_H

/*
 * BATS 线网仿真（只在本进程内存里传）
 *
 * 四节点：源 → 中继1 → 中继2 → 目的。
 * 无预编码、无 socket、无端口。丢包是入队前按概率丢掉。
 *
 * 外码：每个 batch 从 Ψ 抽度数 d，均匀选 d 个源包。
 *       默认是写死的那份。--psi 换一份现成的；--optimize 用本次测到的秩分布自动求一份。
 *       乘 d×M 的随机矩阵 G。G 不进包。
 *       目的端用同一个种子和 batch_id 自己算出 d、下标和 G。
 * 第 j 个发出包的系数向量是 I_M 的第 j 列。payload 已是 (BG) 的第 j 列，
 * 不再乘一次单位阵。
 * 中继：只混合同一个 batch_id 里已经收到的包（没收齐也 recode），
 *       payload 和系数乘同一个随机矩阵。不同 batch 不混。
 * 目的：不 recode。先 BP（rank 等于当前度数就解），解不动再 inactivation。
 *
 * 库（field.c source.c dest.c codec.c）只做编码和译码，见 bats_codec.h。
 * 仿真专用：
 *   link.c    节点之间的内存链路
 *   relay.c   中继。中继1 和中继2 跑同一段 recode
 *   sim.c     源发出包，再串成 源 → 中继1 → 中继2 → 目的
 *   main.c    自检、例子、开销测量、文件还原
 *   psi_opt.c 按论文 (P1) 求 Ψ
 *
 * 用法：
 *   make && ./bats_sim [--optimize | --psi 权重文件] [输入文件] [输出文件]
 *   build.bat              （Windows，用本目录里的 tcc）
 *   --optimize             用本次测量的目的端秩个数按论文 (P1) 求 Ψ，再测一次
 *   psi_opt                单独跑同一套优化，也可以手喂秩个数
 * 还原出的文件与原文一致时返回 0，否则返回 1。
 * 另外对 K=128 从少到多搜索刚好能还原的 batch 数，并打印每跳秩分布和相对 K/h_bar 的开销。
 */

#include "bats_internal.h"

#define LOSS_PERCENT 5
#define FILE_T 32

extern const int PSI_DEG[M];
extern const int PSI_W[M];
extern const int PSI_SUM;

/* 1 表示读入成功。文件第一行是最大度数 D，后面是度数 1..D 的整数权重。 */
int psi_load(const char *path);
/* weights[1..max_degree] 是整数权重，weights[0] 不用。 */
int psi_set(const int *weights, int max_degree);
int psi_is_custom(void);
void psi_print(FILE *fp);
void psi_reset(void);

/* 按 (P1) 从秩个数（秩 0..ncounts-1）求整数权重。
   max_degree<=0 时用论文的度数上限。weights_out 可为 NULL。
   weights_path 可为 NULL。成功返回 1，并把度数上限写入 *D_out。 */
int psi_optimize(const int *counts, int ncounts, int max_degree, double eta_bar, int samples,
                 double q, int scale, int *weights_out, int weights_cap, int *D_out,
                 const char *weights_path);

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

/* 收到的系数向量排成 r×M。整批丢失为 0。 */
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

void eqs_from_packets(BatchEq *eq, const Vec *got, uint64_t code_seed, uint32_t batch_id, int K,
                      int T);

/* 目的：不 recode。把收到的包收成方程，并释放 got。r==0 时不释放（本来就是空的）。 */
void dest_collect(BatchEq *eq, Vec *got, uint64_t code_seed, uint32_t batch, int K, int T,
                  int verbose, int trace, int *rank_out);

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
