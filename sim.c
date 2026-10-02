/* 源 → 中继1 → 中继2 → 目的。每跳之间是一条按概率丢包的内存链路。 */
#include "bats.h"

static int packet_should_print(int verbose, int trace_batch, uint32_t batch_id) {
    if (!verbose) {
        return 0;
    }
    if (trace_batch < 0) {
        return 1;
    }
    return (int)batch_id == trace_batch;
}

SimResult simulate_ex(const uint8_t *src, int K, int T, uint64_t code_seed,
                      uint64_t loss_seed, uint64_t recode_seed, int loss_percent,
                      int n_batches, int verbose, int trace_batch, uint8_t *dst,
                      int *rank_r1, int *rank_r2, int *rank_dst, Vec **dest_keep) {
    SimResult result;
    Rng loss;
    Rng recode;
    BatchEq *eqs = NULL;
    int n_eq = 0;
    int b;
    memset(&result, 0, sizeof(result));
    loss.s = loss_seed;
    recode.s = recode_seed;
    /* 按 batch_id 分桶，避免不同 batch 混在一次 recode 里。 */
    {
        Vec *to_relay1 = xmalloc((size_t)n_batches * sizeof(Vec));
        Vec *to_relay2 = xmalloc((size_t)n_batches * sizeof(Vec));
        Vec *to_dest = xmalloc((size_t)n_batches * sizeof(Vec));
        const char *hop_name[3] = {"[链路 源→中继1]", "[链路 中继1→中继2]", "[链路 中继2→目的]"};
        memset(to_relay1, 0, (size_t)n_batches * sizeof(Vec));
        memset(to_relay2, 0, (size_t)n_batches * sizeof(Vec));
        memset(to_dest, 0, (size_t)n_batches * sizeof(Vec));

        if (verbose) {
            printf("源字节:");
            for (b = 0; b < K; b++) {
                int t;
                printf(" [%d]=", b);
                for (t = 0; t < T; t++) {
                    printf("%02X", src[(size_t)b * (size_t)T + (size_t)t]);
                }
            }
            printf("\n");
        }

        for (b = 0; b < n_batches; b++) {
            int trace = packet_should_print(verbose, trace_batch, (uint32_t)b);
            source_send(src, K, T, code_seed, b, &to_relay1[b], &loss, loss_percent,
                        &result.link[0], verbose, trace, hop_name[0]);
        }

        for (b = 0; b < n_batches; b++) {
            int trace = packet_should_print(verbose, trace_batch, (uint32_t)b);
            relay_forward(&to_relay1[b], &to_relay2[b], &recode, &loss, loss_percent,
                          &result.link[1], b, T, verbose, trace, "中继1", hop_name[1], 1,
                          rank_r1 ? &rank_r1[b] : NULL);
        }
        free(to_relay1);

        for (b = 0; b < n_batches; b++) {
            int trace = packet_should_print(verbose, trace_batch, (uint32_t)b);
            relay_forward(&to_relay2[b], &to_dest[b], &recode, &loss, loss_percent,
                          &result.link[2], b, T, verbose, trace, "中继2", hop_name[2], 0,
                          rank_r2 ? &rank_r2[b] : NULL);
        }
        free(to_relay2);

        if (dest_keep) {
            if (rank_dst) {
                for (b = 0; b < n_batches; b++) {
                    rank_dst[b] = coeff_rank(to_dest[b].p, to_dest[b].n);
                }
            }
            *dest_keep = to_dest;
            return result;
        }

        eqs = xmalloc((size_t)n_batches * sizeof(BatchEq));
        memset(eqs, 0, (size_t)n_batches * sizeof(BatchEq));
        n_eq = n_batches;
        for (b = 0; b < n_batches; b++) {
            int trace = packet_should_print(verbose, trace_batch, (uint32_t)b);
            dest_collect(&eqs[b], &to_dest[b], code_seed, b, K, T, verbose, trace,
                         rank_dst ? &rank_dst[b] : NULL);
        }
        free(to_dest);

        result.ok = decode(eqs, n_eq, K, T, dst, &result.bp_solved, &result.inactivated);
        if (result.ok) {
            result.ok = memcmp(src, dst, (size_t)K * (size_t)T) == 0;
        }
        for (b = 0; b < n_eq; b++) {
            batch_free(&eqs[b]);
        }
        free(eqs);
    }
    return result;
}

SimResult simulate(const uint8_t *src, int K, int T, uint64_t code_seed,
                   uint64_t loss_seed, uint64_t recode_seed, int loss_percent,
                   int n_batches, int verbose, int trace_batch, uint8_t *dst) {
    return simulate_ex(src, K, T, code_seed, loss_seed, recode_seed, loss_percent,
                       n_batches, verbose, trace_batch, dst, NULL, NULL, NULL, NULL);
}

void source_send(const uint8_t *src, int K, int T, uint64_t code_seed, int batch,
                 Vec *next, Rng *loss, int loss_percent, LinkStat *st,
                 int verbose, int trace, const char *link_name) {
    Plan plan = make_plan(code_seed, (uint32_t)batch, K);
    uint8_t *coeff = xmalloc((size_t)M * (size_t)M);
    uint8_t *payload = xmalloc((size_t)M * (size_t)T);
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
    encode_from_plan(src, T, &plan, coeff, payload);
    for (j = 0; j < M; j++) {
        Packet pkt;
        memset(&pkt, 0, sizeof(pkt));
        pkt.batch_id = (uint32_t)batch;
        memcpy(pkt.coeff, coeff + (size_t)j * (size_t)M, (size_t)M);
        pkt.payload = xmalloc((size_t)T);
        memcpy(pkt.payload, payload + (size_t)j * (size_t)T, (size_t)T);
        link_send(pkt, next, loss, loss_percent, st, trace, link_name, j, T);
    }
    free(coeff);
    free(payload);
    if (verbose && !trace) {
        printf("[源] batch=%d 发出 %d 个包（系数是单位阵的列，细节略）\n", batch, M);
    }
    plan_free(&plan);
}

void eqs_from_packets(BatchEq *eq, const Vec *got, uint64_t code_seed, uint32_t batch_id, int K,
                      int T) {
    int r = got->n;
    uint8_t *coeff = NULL;
    uint8_t *payload = NULL;
    int j;
    if (r > 0) {
        coeff = xmalloc((size_t)r * (size_t)M);
        payload = xmalloc((size_t)r * (size_t)T);
        for (j = 0; j < r; j++) {
            memcpy(coeff + (size_t)j * (size_t)M, got->p[j].coeff, (size_t)M);
            memcpy(payload + (size_t)j * (size_t)T, got->p[j].payload, (size_t)T);
        }
    }
    batch_from_received(eq, coeff, payload, r, code_seed, batch_id, K, T);
    free(coeff);
    free(payload);
}

void dest_collect(BatchEq *eq, Vec *got, uint64_t code_seed, uint32_t batch, int K, int T,
                  int verbose, int trace, int *rank_out) {
    int r = got->n;
    int j;
    if (rank_out) {
        *rank_out = coeff_rank(got->p, r);
    }
    if (r == 0) {
        memset(eq, 0, sizeof(*eq));
        eq->batch_id = batch;
        if (verbose) {
            printf("[目的] batch=%u 没有收到包\n", batch);
        }
        return;
    }
    if (verbose) {
        printf("[目的] batch=%u 收到 %d 个，不再乘随机系数\n", batch, r);
    }
    if (trace) {
        for (j = 0; j < r; j++) {
            print_coeff_payload("    [目的] 收到", got->p[j].batch_id, j, got->p[j].coeff,
                                got->p[j].payload, T, "");
        }
    }
    eqs_from_packets(eq, got, code_seed, batch, K, T);
    vec_free(got);
}

int batches_needed(int K) {
    enum { N_HOPS = 3, LN2_MILLI = 693, COVER_SLACK_MILLI = 9000 };
    enum { CUSTOM_BATCH_FACTOR = 6, CUSTOM_BATCH_FLOOR = 32 };
    int log2 = 0;
    int x = K;
    int hop;
    int keep;
    long long deg_sum;
    long long weight_sum;
    long long num;
    long long den;
    long long recv;
    long long percent_scale;
    int n_cover;
    int n_info;
    int n;
    if (K < 1) {
        return 1;
    }
    while (x > 1) {
        x >>= 1;
        log2++;
    }
    weight_sum = psi_weight_sum();
    deg_sum = psi_degree_moment();
    if (deg_sum < 1) {
        deg_sum = 1;
    }
    /* ceil((ln(K)+9) * K / 平均度数)。ln(K) 用 (floor(log2(K))+1)*ln(2) 上界。 */
    num = ((long long)(log2 + 1) * LN2_MILLI + COVER_SLACK_MILLI) * (long long)K * weight_sum;
    den = 1000LL * deg_sum;
    n_cover = (int)((num + den - 1) / den);
    /* 三跳之后每批大约还剩 M*((100-丢包)/100)^3 个包，信息量上界再留一倍。 */
    keep = 100 - LOSS_PERCENT;
    if (keep < 1) {
        keep = 1;
    }
    recv = M;
    percent_scale = 1;
    for (hop = 0; hop < N_HOPS; hop++) {
        recv *= keep;
        percent_scale *= 100;
    }
    num = 2LL * (long long)K * percent_scale;
    n_info = (int)((num + recv - 1) / recv);
    n = n_cover > n_info ? n_cover : n_info;
    if (n < n_info * 2) {
        n = n_info * 2;
    }
    /* 有限 K 比上面的渐近上界费 batch。自定义 Ψ 把搜索上界放到 6K。 */
    if (psi_is_custom()) {
        int wide = CUSTOM_BATCH_FACTOR * K;
        if (wide < CUSTOM_BATCH_FLOOR) {
            wide = CUSTOM_BATCH_FLOOR;
        }
        if (n < wide) {
            n = wide;
        }
    }
    if (n < 1) {
        n = 1;
    }
    return n;
}
