/* 节点之间的内存链路。丢包发生在入队之前。 */
#include "bats.h"

void vec_push(Vec *v, Packet pkt) {
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 8;
        v->p = xrealloc(v->p, (size_t)v->cap * sizeof(Packet));
    }
    v->p[v->n++] = pkt;
}

void vec_free(Vec *v) {
    int i;
    for (i = 0; i < v->n; i++) {
        free(v->p[i].payload);
    }
    free(v->p);
    v->p = NULL;
    v->n = 0;
    v->cap = 0;
}

void print_coeff_payload(const char *tag, uint32_t batch_id, int j,
                                const uint8_t coeff[M], const uint8_t *payload,
                                int T, const char *note) {
    int i;
    printf("%s batch=%u j=%d 系数=[", tag, batch_id, j);
    for (i = 0; i < M; i++) {
        printf("%02X%s", coeff[i], (i + 1 < M) ? "," : "");
    }
    printf("] payload=");
    for (i = 0; i < T; i++) {
        printf("%02X", payload[i]);
    }
    if (note && note[0]) {
        printf(" %s", note);
    }
    printf("\n");
}

int link_send(Packet pkt, Vec *dst, Rng *loss, int loss_percent, LinkStat *st,
                     int do_print, const char *link, int j, int T) {
    int drop = rng_below(loss, 100) < (uint32_t)loss_percent;
    st->sent++;
    if (do_print) {
        print_coeff_payload(link, pkt.batch_id, j, pkt.coeff, pkt.payload, T,
                            drop ? "丢弃" : "送达");
    }
    if (drop) {
        st->dropped++;
        free(pkt.payload);
        return 0;
    }
    st->delivered++;
    vec_push(dst, pkt);
    return 1;
}

int coeff_rank(const Packet *pkts, int r) {
    uint8_t *A;
    int rank;
    int i;
    if (r <= 0) {
        return 0;
    }
    A = xmalloc((size_t)r * (size_t)M);
    for (i = 0; i < r; i++) {
        memcpy(A + (size_t)i * (size_t)M, pkts[i].coeff, (size_t)M);
    }
    rank = coeff_rank_matrix(A, r);
    free(A);
    return rank;
}
