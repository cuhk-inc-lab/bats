/* 源节点：按 Ψ 抽一个 batch，生成 G，发出 M 个包。G 不进包。 */
#include "bats.h"

/* 所有 batch 共用这一份度数分布。权重之和 = 256。
   度数:  1   2   3   4   5   6   7   8
   权重: 12  30  40  44  40  36  30  24
   平均度数 1186/256 ≈ 4.63。最大度数 = M，一个满秩 batch 有机会被 BP 直接解开。 */
const int PSI_DEG[M] = {1, 2, 3, 4, 5, 6, 7, 8};
const int PSI_W[M] = {12, 30, 40, 44, 40, 36, 30, 24};
const int PSI_SUM = 256;

static int sample_degree(Rng *rng) {
    int x = (int)rng_below(rng, (uint32_t)PSI_SUM);
    int acc = 0;
    for (int i = 0; i < M; i++) {
        acc += PSI_W[i];
        if (x < acc) {
            return PSI_DEG[i];
        }
    }
    return PSI_DEG[M - 1];
}

Plan make_plan(uint64_t seed, uint32_t batch_id, int K) {
    Plan plan;
    Rng rng = rng_for_batch(seed, batch_id);
    int d = sample_degree(&rng);
    int *pool;
    int i;
    if (d > K) {
        d = K;
    }
    plan.d = d;
    plan.sel = xmalloc((size_t)d * sizeof(int));
    plan.G = xmalloc((size_t)d * (size_t)M);
    pool = xmalloc((size_t)K * sizeof(int));
    for (i = 0; i < K; i++) {
        pool[i] = i;
    }
    for (i = 0; i < d; i++) {
        int j = i + (int)rng_below(&rng, (uint32_t)(K - i));
        int tmp = pool[i];
        pool[i] = pool[j];
        pool[j] = tmp;
        plan.sel[i] = pool[i];
    }
    for (i = 0; i < d * M; i++) {
        plan.G[i] = rng_byte(&rng);
    }
    free(pool);
    return plan;
}

void plan_free(Plan *plan) {
    free(plan->sel);
    free(plan->G);
    plan->sel = NULL;
    plan->G = NULL;
}

int batches_needed(int K) {
    int log2 = 0;
    int x = K;
    long long num;
    long long den;
    int n_cover;
    int n_info;
    int n;
    while (x > 1) {
        x >>= 1;
        log2++;
    }
    /* (ln(K)+9) * K / d_avg，d_avg = 1186/256。ln 用 (floor(log2)+1)*ln(2) 上界。 */
    num = ((long long)(log2 + 1) * 693 + 9000) * (long long)K * 256;
    den = 1000LL * 1186;
    n_cover = (int)((num + den - 1) / den);
    /* 目的端每批大约收到 M*(0.95^3) 个包，再留一倍余量。 */
    n_info = (2 * K * 100 + 679) / 680;
    n = n_cover > n_info ? n_cover : n_info;
    if (n < 1) {
        n = 1;
    }
    return n;
}

void source_send(const uint8_t *src, int K, int T, uint64_t code_seed, int batch,
                 Vec *next, Rng *loss, int loss_percent, LinkStat *st,
                 int verbose, int trace, const char *link_name) {
    Plan plan = make_plan(code_seed, (uint32_t)batch, K);
    int j;
    if (trace) {
        int a;
        int m;
        printf("[源] batch=%d 只是箱子编号  d=%d 选中源包", batch, plan.d);
        for (a = 0; a < plan.d; a++) {
            printf(" %d", plan.sel[a]);
        }
        printf("\n[源] batch=%d 的 G（%d×%d）不放进包，目的端用同一种子复原\n",
               batch, plan.d, M);
        for (a = 0; a < plan.d; a++) {
            printf("    对应源包 %d:", plan.sel[a]);
            for (m = 0; m < M; m++) {
                printf(" %02X", plan.G[(size_t)a * M + (size_t)m]);
            }
            printf("\n");
        }
    }
    for (j = 0; j < M; j++) {
        Packet pkt;
        int t;
        int a;
        memset(&pkt, 0, sizeof(pkt));
        pkt.batch_id = (uint32_t)batch;
        pkt.coeff[j] = 1; /* I_M 的第 j 列 */
        pkt.payload = xmalloc((size_t)T);
        memset(pkt.payload, 0, (size_t)T);
        for (a = 0; a < plan.d; a++) {
            uint8_t g = plan.G[(size_t)a * M + (size_t)j];
            const uint8_t *sp;
            if (!g) {
                continue;
            }
            sp = src + (size_t)plan.sel[a] * (size_t)T;
            for (t = 0; t < T; t++) {
                pkt.payload[t] ^= gf_mul(g, sp[t]);
            }
        }
        link_send(pkt, next, loss, loss_percent, st, trace, link_name, j, T);
    }
    if (verbose && !trace) {
        printf("[源] batch=%d 发出 %d 个包（系数是单位阵的列，细节略）\n", batch, M);
    }
    plan_free(&plan);
}
