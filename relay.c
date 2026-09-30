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
    for (j = 0; j < M; j++) {
        Packet pkt;
        int nz = 0;
        uint8_t *phi_col = xmalloc((size_t)r);
        memset(&pkt, 0, sizeof(pkt));
        do {
            nz = 0;
            for (s = 0; s < r; s++) {
                phi_col[s] = rng_byte(recode);
                if (phi_col[s]) {
                    nz = 1;
                }
            }
        } while (!nz);
        pkt.batch_id = (uint32_t)batch;
        pkt.payload = xmalloc((size_t)T);
        memset(pkt.payload, 0, (size_t)T);
        for (s = 0; s < r; s++) {
            uint8_t a = phi_col[s];
            int m;
            int t;
            if (!a) {
                continue;
            }
            for (m = 0; m < M; m++) {
                pkt.coeff[m] ^= gf_mul(a, in->p[s].coeff[m]);
            }
            for (t = 0; t < T; t++) {
                pkt.payload[t] ^= gf_mul(a, in->p[s].payload[t]);
            }
        }
        free(phi_col);
        if (trace) {
            print_coeff_payload(send_tag, pkt.batch_id, j, pkt.coeff, pkt.payload, T, "");
        }
        link_send(pkt, next, loss, loss_percent, st, trace, link_name, j, T);
    }
    vec_free(in);
}
