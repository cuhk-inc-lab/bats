/* 编码、recode、译码。仿真和外部程序都走这里。 */
#include "bats_internal.h"

void bats_init(void) {
    gf_init();
}

void encode_from_plan(const uint8_t *src, int T, const Plan *plan, uint8_t *coeff, uint8_t *payload) {
    int j;
    for (j = 0; j < M; j++) {
        int a;
        uint8_t *cj = coeff + (size_t)j * (size_t)M;
        uint8_t *pj = payload + (size_t)j * (size_t)T;
        memset(cj, 0, (size_t)M);
        cj[j] = 1; /* I_M 的第 j 列 */
        memset(pj, 0, (size_t)T);
        for (a = 0; a < plan->d; a++) {
            uint8_t g = plan->G[(size_t)a * (size_t)M + (size_t)j];
            const uint8_t *sp;
            if (!g) {
                continue;
            }
            sp = src + (size_t)plan->sel[a] * (size_t)T;
            gf_axpy(pj, sp, g, (size_t)T);
        }
    }
}

int bats_encode(const uint8_t *src, int K, int T, uint64_t code_seed, uint32_t batch_id,
                uint8_t *coeff, uint8_t *payload) {
    Plan plan;
    if (K < 1 || T < 1 || src == NULL || coeff == NULL || payload == NULL) {
        return -1;
    }
    plan = make_plan(code_seed, batch_id, K);
    encode_from_plan(src, T, &plan, coeff, payload);
    plan_free(&plan);
    return 0;
}

int bats_recode(const uint8_t *coeff_in, const uint8_t *payload_in, int n_in, int T,
                uint64_t *recode_state, uint8_t *coeff_out, uint8_t *payload_out,
                int *n_out) {
    Rng recode;
    uint8_t coeff_stack[M * M];
    uint8_t payload_stack[M * 2048];
    uint8_t phi_stack[128];
    uint8_t *coeff_tmp;
    uint8_t *payload_tmp;
    uint8_t *phi_col;
    int use_heap;
    int j;
    if (n_out) {
        *n_out = 0;
    }
    if (n_in < 0 || T < 1 || recode_state == NULL) {
        return -1;
    }
    if (n_in == 0) {
        return 0;
    }
    if (coeff_in == NULL || payload_in == NULL || coeff_out == NULL || payload_out == NULL) {
        return -1;
    }
    recode.s = *recode_state;
    use_heap = (T > 2048 || n_in > 128);
    if (use_heap) {
        coeff_tmp = xmalloc((size_t)M * (size_t)M);
        payload_tmp = xmalloc((size_t)M * (size_t)T);
        phi_col = xmalloc((size_t)n_in);
    } else {
        coeff_tmp = coeff_stack;
        payload_tmp = payload_stack;
        phi_col = phi_stack;
    }
    for (j = 0; j < M; j++) {
        int nz = 0;
        int s;
        uint8_t *cj = coeff_tmp + (size_t)j * (size_t)M;
        uint8_t *pj = payload_tmp + (size_t)j * (size_t)T;
        memset(cj, 0, (size_t)M);
        memset(pj, 0, (size_t)T);
        do {
            nz = 0;
            for (s = 0; s < n_in; s++) {
                phi_col[s] = rng_byte(&recode);
                if (phi_col[s]) {
                    nz = 1;
                }
            }
        } while (!nz);
        for (s = 0; s < n_in; s++) {
            uint8_t a = phi_col[s];
            int m;
            const uint8_t *cs;
            const uint8_t *ps;
            if (!a) {
                continue;
            }
            cs = coeff_in + (size_t)s * (size_t)M;
            ps = payload_in + (size_t)s * (size_t)T;
            for (m = 0; m < M; m++) {
                cj[m] ^= gf_mul(a, cs[m]);
            }
            gf_axpy(pj, ps, a, (size_t)T);
        }
    }
    memcpy(coeff_out, coeff_tmp, (size_t)M * (size_t)M);
    memcpy(payload_out, payload_tmp, (size_t)M * (size_t)T);
    if (use_heap) {
        free(phi_col);
        free(coeff_tmp);
        free(payload_tmp);
    }
    *recode_state = recode.s;
    if (n_out) {
        *n_out = M;
    }
    return 0;
}

int bats_decode(const BatsSymbol *symbols, int n_symbols, int K, int T, uint64_t code_seed,
                uint8_t *dst, int *n_bp, int *n_inact) {
    uint32_t *ids = NULL;
    int *counts = NULL;
    int n_batches = 0;
    BatchEq *eqs = NULL;
    int i;
    int ok;
    if (n_bp) {
        *n_bp = 0;
    }
    if (n_inact) {
        *n_inact = 0;
    }
    if (K < 1 || T < 1 || dst == NULL || n_symbols < 0) {
        return -1;
    }
    if (n_symbols == 0) {
        return 0;
    }
    if (symbols == NULL) {
        return -1;
    }
    ids = xmalloc((size_t)n_symbols * sizeof(uint32_t));
    counts = xmalloc((size_t)n_symbols * sizeof(int));
    for (i = 0; i < n_symbols; i++) {
        int b;
        int pos = -1;
        if (symbols[i].payload == NULL) {
            free(ids);
            free(counts);
            return -1;
        }
        for (b = 0; b < n_batches; b++) {
            if (ids[b] == symbols[i].batch_id) {
                pos = b;
                break;
            }
        }
        if (pos < 0) {
            pos = n_batches;
            ids[n_batches] = symbols[i].batch_id;
            counts[n_batches] = 0;
            n_batches++;
        }
        counts[pos]++;
    }
    for (i = 1; i < n_batches; i++) {
        uint32_t id = ids[i];
        int count = counts[i];
        int j = i;
        while (j > 0 && ids[j - 1] > id) {
            ids[j] = ids[j - 1];
            counts[j] = counts[j - 1];
            j--;
        }
        ids[j] = id;
        counts[j] = count;
    }
    eqs = xmalloc((size_t)n_batches * sizeof(BatchEq));
    memset(eqs, 0, (size_t)n_batches * sizeof(BatchEq));
    for (i = 0; i < n_batches; i++) {
        uint8_t *coeff = xmalloc((size_t)counts[i] * (size_t)M);
        uint8_t *payload = xmalloc((size_t)counts[i] * (size_t)T);
        int filled = 0;
        int s;
        for (s = 0; s < n_symbols; s++) {
            if (symbols[s].batch_id != ids[i]) {
                continue;
            }
            memcpy(coeff + (size_t)filled * (size_t)M, symbols[s].coeff, (size_t)M);
            memcpy(payload + (size_t)filled * (size_t)T, symbols[s].payload, (size_t)T);
            filled++;
        }
        batch_from_received(&eqs[i], coeff, payload, filled, code_seed, ids[i], K, T);
        free(coeff);
        free(payload);
    }
    ok = decode(eqs, n_batches, K, T, dst, n_bp, n_inact);
    for (i = 0; i < n_batches; i++) {
        batch_free(&eqs[i]);
    }
    free(eqs);
    free(ids);
    free(counts);
    return ok;
}
