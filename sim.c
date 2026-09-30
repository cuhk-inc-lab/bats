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
