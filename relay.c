/* 中继节点：只混合同一个 batch 里已经收到的包。中继1 和中继2 跑这一段。 */
#include "bats.h"

void relay_forward(Vec *in, Vec *next, Rng *recode, Rng *loss, int loss_percent,
                   LinkStat *st, int batch, int T, int verbose, int trace,
                   const char *name, const char *link_name, int mention_m,
                   int *rank_out) {
    int r = in->n;
    int j;
    int s;
    char recv_tag[64];
    char send_tag[64];
    if (rank_out) {
        *rank_out = coeff_rank(in->p, r);
    }
    if (r == 0) {
        if (verbose) {
            printf("[%s] batch=%d 一个包都没收到，不 recode\n", name, batch);
        }
        return;
    }
    if (verbose) {
        if (mention_m) {
            printf("[%s] batch=%d 收到 %d 个（M=%d，未收齐也 recode）\n", name, batch, r, M);
        } else {
            printf("[%s] batch=%d 收到 %d 个（未收齐也 recode）\n", name, batch, r);
        }
    }
    snprintf(recv_tag, sizeof(recv_tag), "    [%s] 收到", name);
    snprintf(send_tag, sizeof(send_tag), "    [%s] recode 发出", name);
    if (trace) {
        for (s = 0; s < r; s++) {
            print_coeff_payload(recv_tag, in->p[s].batch_id, s, in->p[s].coeff,
                                in->p[s].payload, T, "");
        }
    }
    {
        uint8_t *coeff_in = xmalloc((size_t)r * (size_t)M);
        uint8_t *payload_in = xmalloc((size_t)r * (size_t)T);
        uint8_t *coeff_out = xmalloc((size_t)M * (size_t)M);
        uint8_t *payload_out = xmalloc((size_t)M * (size_t)T);
        int n_out = 0;
        for (s = 0; s < r; s++) {
            memcpy(coeff_in + (size_t)s * (size_t)M, in->p[s].coeff, (size_t)M);
            memcpy(payload_in + (size_t)s * (size_t)T, in->p[s].payload, (size_t)T);
        }
        bats_recode(coeff_in, payload_in, r, T, &recode->s, coeff_out, payload_out, &n_out);
        for (j = 0; j < n_out; j++) {
            Packet pkt;
            memset(&pkt, 0, sizeof(pkt));
            pkt.batch_id = (uint32_t)batch;
            memcpy(pkt.coeff, coeff_out + (size_t)j * (size_t)M, (size_t)M);
            pkt.payload = xmalloc((size_t)T);
            memcpy(pkt.payload, payload_out + (size_t)j * (size_t)T, (size_t)T);
            if (trace) {
                print_coeff_payload(send_tag, pkt.batch_id, j, pkt.coeff, pkt.payload, T, "");
            }
            link_send(pkt, next, loss, loss_percent, st, trace, link_name, j, T);
        }
        free(coeff_in);
        free(payload_in);
        free(coeff_out);
        free(payload_out);
    }
    vec_free(in);
}
