/* 目的节点：不 recode。先 BP，解不动再 inactivation。 */
#include "bats.h"

enum { VAR_UNKNOWN = 0, VAR_INACT = 1, VAR_SOLVED = 2 };

typedef struct {
    int status;
    int inact_id;
    uint8_t *cons;
    int *dep_id;
    uint8_t *dep_c;
    int dep_n;
    int dep_cap;
} Var;

static void batch_push_row(BatchEq *b, int var, const uint8_t *row) {
    int all0 = 1;
    int j;
    for (j = 0; j < b->n_recv; j++) {
        if (row[j]) {
            all0 = 0;
            break;
        }
    }
    if (all0) {
        return;
    }
    if (b->n_nb == b->cap_nb) {
        b->cap_nb = b->cap_nb ? b->cap_nb * 2 : 4;
        b->nb = xrealloc(b->nb, (size_t)b->cap_nb * sizeof(int));
        b->gamma = xrealloc(b->gamma, (size_t)b->cap_nb * (size_t)b->n_recv);
    }
    b->nb[b->n_nb] = var;
    memcpy(b->gamma + (size_t)b->n_nb * (size_t)b->n_recv, row, (size_t)b->n_recv);
    b->n_nb++;
}

static void batch_remove_row(BatchEq *b, int k) {
    int nR = b->n_recv;
    int tail = b->n_nb - k - 1;
    if (tail > 0) {
        memmove(b->nb + k, b->nb + k + 1, (size_t)tail * sizeof(int));
        memmove(b->gamma + (size_t)k * (size_t)nR,
                b->gamma + (size_t)(k + 1) * (size_t)nR,
                (size_t)tail * (size_t)nR);
    }
    b->n_nb--;
}

static int batch_find(const BatchEq *b, int var) {
    int k;
    for (k = 0; k < b->n_nb; k++) {
        if (b->nb[k] == var) {
            return k;
        }
    }
    return -1;
}

static void batch_accumulate(BatchEq *b, int var, const uint8_t *delta) {
    int k = batch_find(b, var);
    int j;
    int all0 = 1;
    if (k >= 0) {
        for (j = 0; j < b->n_recv; j++) {
            b->gamma[(size_t)k * (size_t)b->n_recv + (size_t)j] ^= delta[j];
        }
        for (j = 0; j < b->n_recv; j++) {
            if (b->gamma[(size_t)k * (size_t)b->n_recv + (size_t)j]) {
                all0 = 0;
                break;
            }
        }
        if (all0) {
            batch_remove_row(b, k);
        }
        return;
    }
    batch_push_row(b, var, delta);
}

static void strip_zero_rows(BatchEq *b) {
    int nR = b->n_recv;
    int w = 0;
    int k;
    for (k = 0; k < b->n_nb; k++) {
        int nz = 0;
        int j;
        for (j = 0; j < nR; j++) {
            if (b->gamma[(size_t)k * (size_t)nR + (size_t)j]) {
                nz = 1;
                break;
            }
        }
        if (!nz) {
            continue;
        }
        if (w != k) {
            b->nb[w] = b->nb[k];
            memcpy(b->gamma + (size_t)w * (size_t)nR,
                   b->gamma + (size_t)k * (size_t)nR, (size_t)nR);
        }
        w++;
    }
    b->n_nb = w;
}

void batch_free(BatchEq *b) {
    free(b->nb);
    free(b->gamma);
    free(b->Y);
    b->nb = NULL;
    b->gamma = NULL;
    b->Y = NULL;
}

static void var_clear_deps(Var *v) {
    free(v->dep_id);
    free(v->dep_c);
    v->dep_id = NULL;
    v->dep_c = NULL;
    v->dep_n = 0;
    v->dep_cap = 0;
}

static void var_add_dep(Var *v, int id, uint8_t c) {
    int i;
    if (c == 0) {
        return;
    }
    for (i = 0; i < v->dep_n; i++) {
        if (v->dep_id[i] == id) {
            v->dep_c[i] ^= c;
            if (v->dep_c[i] == 0) {
                v->dep_n--;
                v->dep_id[i] = v->dep_id[v->dep_n];
                v->dep_c[i] = v->dep_c[v->dep_n];
            }
            return;
        }
    }
    if (v->dep_n == v->dep_cap) {
        v->dep_cap = v->dep_cap ? v->dep_cap * 2 : 4;
        v->dep_id = xrealloc(v->dep_id, (size_t)v->dep_cap * sizeof(int));
        v->dep_c = xrealloc(v->dep_c, (size_t)v->dep_cap);
    }
    v->dep_id[v->dep_n] = id;
    v->dep_c[v->dep_n] = c;
    v->dep_n++;
}

static void eliminate_var(BatchEq *b, int var, Var *v, int T) {
    int k = batch_find(b, var);
    int nR;
    uint8_t *row;
    int j;
    int di;
    if (k < 0) {
        return;
    }
    nR = b->n_recv;
    row = xmalloc((size_t)nR);
    memcpy(row, b->gamma + (size_t)k * (size_t)nR, (size_t)nR);
    batch_remove_row(b, k);
    for (j = 0; j < nR; j++) {
        int t;
        if (!row[j]) {
            continue;
        }
        for (t = 0; t < T; t++) {
            b->Y[(size_t)j * (size_t)T + (size_t)t] ^= gf_mul(row[j], v->cons[t]);
        }
    }
    for (di = 0; di < v->dep_n; di++) {
        uint8_t *delta = xmalloc((size_t)nR);
        for (j = 0; j < nR; j++) {
            delta[j] = gf_mul(row[j], v->dep_c[di]);
        }
        batch_accumulate(b, v->dep_id[di], delta);
        free(delta);
    }
    free(row);
}

static int try_solve_batch(BatchEq *b, Var *vars, int T, BatchEq *all, int n_batches) {
    int nR = b->n_recv;
    int *U;
    int nU = 0;
    int k;
    uint8_t *A;
    uint8_t *Inv;
    int *piv;
    int rank;
    int u;
    int solved_n = 0;
    Var *saved;
    int *solved_var;
    if (nR == 0) {
        return 0;
    }
    strip_zero_rows(b);
    U = xmalloc((size_t)(b->n_nb > 0 ? b->n_nb : 1) * sizeof(int));
    for (k = 0; k < b->n_nb; k++) {
        if (vars[b->nb[k]].status == VAR_UNKNOWN) {
            U[nU++] = k;
        }
    }
    if (nU == 0) {
        free(U);
        return 0;
    }
    A = xmalloc((size_t)nU * (size_t)nR);
    for (u = 0; u < nU; u++) {
        memcpy(A + (size_t)u * (size_t)nR,
               b->gamma + (size_t)U[u] * (size_t)nR, (size_t)nR);
    }
    Inv = xmalloc((size_t)nU * (size_t)nU);
    piv = xmalloc((size_t)nU * sizeof(int));
    rank = row_reduce(A, nU, nR, piv, Inv);
    if (rank != nU) {
        free(U);
        free(A);
        free(Inv);
        free(piv);
        return 0;
    }
    saved = xmalloc((size_t)nU * sizeof(Var));
    solved_var = xmalloc((size_t)nU * sizeof(int));
    memset(saved, 0, (size_t)nU * sizeof(Var));
    for (u = 0; u < nU; u++) {
        int i;
        Var *sv = &saved[u];
        solved_var[u] = b->nb[U[u]];
        sv->cons = xmalloc((size_t)T);
        memset(sv->cons, 0, (size_t)T);
        for (i = 0; i < nU; i++) {
            uint8_t coef = Inv[(size_t)i * (size_t)nU + (size_t)u];
            int col = piv[i];
            int t;
            if (!coef) {
                continue;
            }
            for (t = 0; t < T; t++) {
                sv->cons[t] ^= gf_mul(coef, b->Y[(size_t)col * (size_t)T + (size_t)t]);
            }
        }
        for (k = 0; k < b->n_nb; k++) {
            int is_u = 0;
            uint8_t alpha = 0;
            int uu;
            for (uu = 0; uu < nU; uu++) {
                if (U[uu] == k) {
                    is_u = 1;
                    break;
                }
            }
            if (is_u) {
                continue;
            }
            if (vars[b->nb[k]].status != VAR_INACT) {
                continue;
            }
            for (i = 0; i < nU; i++) {
                int col = piv[i];
                uint8_t coef = Inv[(size_t)i * (size_t)nU + (size_t)u];
                alpha ^= gf_mul(coef, b->gamma[(size_t)k * (size_t)nR + (size_t)col]);
            }
            var_add_dep(sv, vars[b->nb[k]].inact_id, alpha);
        }
        solved_n++;
    }
    for (u = 0; u < solved_n; u++) {
        Var *dstv = &vars[solved_var[u]];
        free(dstv->cons);
        var_clear_deps(dstv);
        dstv->status = VAR_SOLVED;
        dstv->cons = saved[u].cons;
        dstv->dep_id = saved[u].dep_id;
        dstv->dep_c = saved[u].dep_c;
        dstv->dep_n = saved[u].dep_n;
        dstv->dep_cap = saved[u].dep_cap;
        saved[u].cons = NULL;
        saved[u].dep_id = NULL;
        saved[u].dep_c = NULL;
    }
    for (u = 0; u < solved_n; u++) {
        int bi;
        for (bi = 0; bi < n_batches; bi++) {
            eliminate_var(&all[bi], solved_var[u], &vars[solved_var[u]], T);
        }
    }
    free(saved);
    free(solved_var);
    free(U);
    free(A);
    free(Inv);
    free(piv);
    return 1;
}

static int count_status(const Var *vars, int K, int status) {
    int i;
    int n = 0;
    for (i = 0; i < K; i++) {
        if (vars[i].status == status) {
            n++;
        }
    }
    return n;
}

static int choose_inactivate(BatchEq *batches, int n_batches, Var *vars, int K) {
    int *occ = xmalloc((size_t)K * sizeof(int));
    int b;
    int best_b = -1;
    int best_gap = 0x7fffffff;
    int best_var = -1;
    int best_occ = -1;
    memset(occ, 0, (size_t)K * sizeof(int));
    for (b = 0; b < n_batches; b++) {
        int k;
        strip_zero_rows(&batches[b]);
        for (k = 0; k < batches[b].n_nb; k++) {
            int v = batches[b].nb[k];
            if (vars[v].status == VAR_UNKNOWN) {
                occ[v]++;
            }
        }
    }
    for (b = 0; b < n_batches; b++) {
        BatchEq *be = &batches[b];
        int nR = be->n_recv;
        int nU = 0;
        int k;
        uint8_t *A;
        int rank;
        int gap;
        for (k = 0; k < be->n_nb; k++) {
            if (vars[be->nb[k]].status == VAR_UNKNOWN) {
                nU++;
            }
        }
        if (nU == 0) {
            continue;
        }
        A = xmalloc((size_t)nU * (size_t)nR);
        {
            int u = 0;
            for (k = 0; k < be->n_nb; k++) {
                if (vars[be->nb[k]].status == VAR_UNKNOWN) {
                    memcpy(A + (size_t)u * (size_t)nR,
                           be->gamma + (size_t)k * (size_t)nR, (size_t)nR);
                    u++;
                }
            }
        }
        rank = row_reduce(A, nU, nR, NULL, NULL);
        free(A);
        gap = nU - rank;
        if (gap > 0 && gap < best_gap) {
            best_gap = gap;
            best_b = b;
        }
    }
    if (best_b >= 0) {
        BatchEq *be = &batches[best_b];
        int k;
        for (k = 0; k < be->n_nb; k++) {
            int v = be->nb[k];
            if (vars[v].status == VAR_UNKNOWN && occ[v] > best_occ) {
                best_occ = occ[v];
                best_var = v;
            }
        }
    }
    if (best_var < 0) {
        int v;
        for (v = 0; v < K; v++) {
            if (vars[v].status == VAR_UNKNOWN && occ[v] > best_occ) {
                best_occ = occ[v];
                best_var = v;
            }
        }
    }
    free(occ);
    return best_var;
}

typedef struct {
    uint8_t *A;
    uint8_t *rhs;
    int n;
    int cap;
    int nvar;
    int T;
} Sys;

static void sys_push(Sys *s, const uint8_t *coeff, const uint8_t *rhs) {
    int all0 = 1;
    int j;
    int t;
    for (j = 0; j < s->nvar; j++) {
        if (coeff[j]) {
            all0 = 0;
            break;
        }
    }
    if (all0) {
        for (t = 0; t < s->T; t++) {
            if (rhs[t]) {
                s->n = -1; /* 0 = 非零，矛盾 */
                return;
            }
        }
        return;
    }
    if (s->n < 0) {
        return;
    }
    if (s->n == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 8;
        s->A = xrealloc(s->A, (size_t)s->cap * (size_t)s->nvar);
        s->rhs = xrealloc(s->rhs, (size_t)s->cap * (size_t)s->T);
    }
    memcpy(s->A + (size_t)s->n * (size_t)s->nvar, coeff, (size_t)s->nvar);
    memcpy(s->rhs + (size_t)s->n * (size_t)s->T, rhs, (size_t)s->T);
    s->n++;
}

/* 解 A x = rhs。成功则 sol 按变量给出 T 字节。 */
static int solve_dense(uint8_t *A, uint8_t *rhs, int n_eq, int nvar, int T, uint8_t *sol) {
    int *where = xmalloc((size_t)nvar * sizeof(int));
    int row = 0;
    int col;
    int i;
    for (col = 0; col < nvar; col++) {
        where[col] = -1;
    }
    for (col = 0; col < nvar && row < n_eq; col++) {
        int piv = -1;
        uint8_t invp;
        for (i = row; i < n_eq; i++) {
            if (A[(size_t)i * (size_t)nvar + (size_t)col]) {
                piv = i;
                break;
            }
        }
        if (piv < 0) {
            continue;
        }
        if (piv != row) {
            int c;
            for (c = 0; c < nvar; c++) {
                uint8_t tmp = A[(size_t)row * (size_t)nvar + (size_t)c];
                A[(size_t)row * (size_t)nvar + (size_t)c] =
                    A[(size_t)piv * (size_t)nvar + (size_t)c];
                A[(size_t)piv * (size_t)nvar + (size_t)c] = tmp;
            }
            for (c = 0; c < T; c++) {
                uint8_t tmp = rhs[(size_t)row * (size_t)T + (size_t)c];
                rhs[(size_t)row * (size_t)T + (size_t)c] =
                    rhs[(size_t)piv * (size_t)T + (size_t)c];
                rhs[(size_t)piv * (size_t)T + (size_t)c] = tmp;
            }
        }
        invp = gf_inv(A[(size_t)row * (size_t)nvar + (size_t)col]);
        {
            int c;
            for (c = 0; c < nvar; c++) {
                A[(size_t)row * (size_t)nvar + (size_t)c] =
                    gf_mul(A[(size_t)row * (size_t)nvar + (size_t)c], invp);
            }
            for (c = 0; c < T; c++) {
                rhs[(size_t)row * (size_t)T + (size_t)c] =
                    gf_mul(rhs[(size_t)row * (size_t)T + (size_t)c], invp);
            }
        }
        for (i = 0; i < n_eq; i++) {
            uint8_t f;
            int c;
            if (i == row) {
                continue;
            }
            f = A[(size_t)i * (size_t)nvar + (size_t)col];
            if (!f) {
                continue;
            }
            for (c = 0; c < nvar; c++) {
                A[(size_t)i * (size_t)nvar + (size_t)c] ^=
                    gf_mul(f, A[(size_t)row * (size_t)nvar + (size_t)c]);
            }
            for (c = 0; c < T; c++) {
                rhs[(size_t)i * (size_t)T + (size_t)c] ^=
                    gf_mul(f, rhs[(size_t)row * (size_t)T + (size_t)c]);
            }
        }
        where[col] = row;
        row++;
    }
    for (i = 0; i < n_eq; i++) {
        int all0 = 1;
        int t;
        for (col = 0; col < nvar; col++) {
            if (A[(size_t)i * (size_t)nvar + (size_t)col]) {
                all0 = 0;
                break;
            }
        }
        if (!all0) {
            continue;
        }
        for (t = 0; t < T; t++) {
            if (rhs[(size_t)i * (size_t)T + (size_t)t]) {
                free(where);
                return 0;
            }
        }
    }
    for (col = 0; col < nvar; col++) {
        if (where[col] < 0) {
            free(where);
            return 0;
        }
        memcpy(sol + (size_t)col * (size_t)T,
               rhs + (size_t)where[col] * (size_t)T, (size_t)T);
    }
    free(where);
    return 1;
}

int decode(BatchEq *batches, int n_batches, int K, int T, uint8_t *dst,
                  int *n_bp_out, int *n_inact_out) {
    Var *vars = xmalloc((size_t)K * sizeof(Var));
    uint8_t *inact_val = NULL;
    int n_inact = 0;
    int n_bp = 0;
    int ok = 0;
    int i;
    memset(vars, 0, (size_t)K * sizeof(Var));
    memset(dst, 0, (size_t)K * (size_t)T);
    for (;;) {
        int progressed = 0;
        int b;
        int unknown;
        for (b = 0; b < n_batches; b++) {
            if (try_solve_batch(&batches[b], vars, T, batches, n_batches)) {
                progressed = 1;
                break;
            }
        }
        if (progressed) {
            continue;
        }
        unknown = count_status(vars, K, VAR_UNKNOWN);
        if (unknown == 0) {
            break;
        }
        {
            int v = choose_inactivate(batches, n_batches, vars, K);
            if (v < 0) {
                goto done;
            }
            vars[v].status = VAR_INACT;
            vars[v].inact_id = n_inact;
            n_inact++;
        }
    }
    n_bp = count_status(vars, K, VAR_SOLVED);
    if (n_inact > 0) {
        Sys sys;
        uint8_t *coeff;
        int b;
        memset(&sys, 0, sizeof(sys));
        sys.nvar = n_inact;
        sys.T = T;
        coeff = xmalloc((size_t)n_inact);
        for (b = 0; b < n_batches; b++) {
            BatchEq *be = &batches[b];
            int j;
            int guard;
            for (guard = 0; guard < K + 2; guard++) {
                int k;
                int found = 0;
                strip_zero_rows(be);
                for (k = 0; k < be->n_nb; k++) {
                    if (vars[be->nb[k]].status == VAR_SOLVED) {
                        eliminate_var(be, be->nb[k], &vars[be->nb[k]], T);
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    break;
                }
            }
            for (j = 0; j < be->n_recv; j++) {
                int k;
                memset(coeff, 0, (size_t)n_inact);
                for (k = 0; k < be->n_nb; k++) {
                    Var *v = &vars[be->nb[k]];
                    if (v->status != VAR_INACT) {
                        free(coeff);
                        goto done;
                    }
                    coeff[v->inact_id] ^= be->gamma[(size_t)k * (size_t)be->n_recv + (size_t)j];
                }
                sys_push(&sys, coeff, be->Y + (size_t)j * (size_t)T);
                if (sys.n < 0) {
                    free(coeff);
                    free(sys.A);
                    free(sys.rhs);
                    goto done;
                }
            }
        }
        free(coeff);
        inact_val = xmalloc((size_t)n_inact * (size_t)T);
        if (sys.n == 0 || !solve_dense(sys.A, sys.rhs, sys.n, n_inact, T, inact_val)) {
            free(sys.A);
            free(sys.rhs);
            goto done;
        }
        free(sys.A);
        free(sys.rhs);
    }
    for (i = 0; i < K; i++) {
        int t;
        int di;
        if (vars[i].status == VAR_INACT) {
            memcpy(dst + (size_t)i * (size_t)T,
                   inact_val + (size_t)vars[i].inact_id * (size_t)T, (size_t)T);
        } else if (vars[i].status == VAR_SOLVED) {
            for (di = 0; di < vars[i].dep_n; di++) {
                int id = vars[i].dep_id[di];
                uint8_t c = vars[i].dep_c[di];
                for (t = 0; t < T; t++) {
                    vars[i].cons[t] ^= gf_mul(c, inact_val[(size_t)id * (size_t)T + (size_t)t]);
                }
            }
            memcpy(dst + (size_t)i * (size_t)T, vars[i].cons, (size_t)T);
        } else {
            goto done;
        }
    }
    ok = 1;
done:
    if (n_bp_out) {
        *n_bp_out = n_bp;
    }
    if (n_inact_out) {
        *n_inact_out = n_inact;
    }
    for (i = 0; i < K; i++) {
        free(vars[i].cons);
        var_clear_deps(&vars[i]);
    }
    free(vars);
    free(inact_val);
    return ok;
}

int inactivation_self_test(void) {
    /* 两个度数为 2、秩为 1 的 batch。BP 解不动，必须 inactivation。 */
    const uint8_t x = 4;
    const uint8_t y = 6;
    uint8_t y0 = (uint8_t)(x ^ y);
    uint8_t y1 = (uint8_t)(x ^ gf_mul(3, y));
    BatchEq batches[2];
    uint8_t row[1];
    uint8_t dst[2];
    int n_bp = 0;
    int n_inact = 0;
    int ok;
    memset(batches, 0, sizeof(batches));
    batches[0].n_recv = 1;
    batches[0].Y = xmalloc(1);
    batches[0].Y[0] = y0;
    batches[0].batch_id = 0;
    row[0] = 1;
    batch_push_row(&batches[0], 0, row);
    row[0] = 1;
    batch_push_row(&batches[0], 1, row);
    batches[1].n_recv = 1;
    batches[1].Y = xmalloc(1);
    batches[1].Y[0] = y1;
    batches[1].batch_id = 1;
    row[0] = 1;
    batch_push_row(&batches[1], 0, row);
    row[0] = 3;
    batch_push_row(&batches[1], 1, row);
    ok = decode(batches, 2, 2, 1, dst, &n_bp, &n_inact);
    batch_free(&batches[0]);
    batch_free(&batches[1]);
    if (!ok || dst[0] != x || dst[1] != y || n_inact < 1) {
        fprintf(stderr, "inactivation 自检失败 ok=%d dst=%02X %02X inact=%d\n",
                ok, dst[0], dst[1], n_inact);
        return 0;
    }
    return 1;
}

void eqs_from_packets(BatchEq *eq, const Vec *got, uint64_t code_seed, int batch_id, int K, int T) {
    int r = got->n;
    Plan plan;
    int a;
    int j;
    memset(eq, 0, sizeof(*eq));
    eq->batch_id = batch_id;
    if (r == 0) {
        return;
    }
    plan = make_plan(code_seed, (uint32_t)batch_id, K);
    eq->n_recv = r;
    eq->Y = xmalloc((size_t)r * (size_t)T);
    for (j = 0; j < r; j++) {
        memcpy(eq->Y + (size_t)j * (size_t)T, got->p[j].payload, (size_t)T);
    }
    for (a = 0; a < plan.d; a++) {
        uint8_t *row = xmalloc((size_t)r);
        for (j = 0; j < r; j++) {
            uint8_t acc = 0;
            int m;
            for (m = 0; m < M; m++) {
                acc ^= gf_mul(plan.G[(size_t)a * (size_t)M + (size_t)m], got->p[j].coeff[m]);
            }
            row[j] = acc;
        }
        batch_push_row(eq, plan.sel[a], row);
        free(row);
    }
    plan_free(&plan);
}

void dest_collect(BatchEq *eq, Vec *got, uint64_t code_seed, int batch, int K, int T,
                  int verbose, int trace, int *rank_out) {
    int r = got->n;
    Plan plan;
    int a;
    int j;
    eq->batch_id = batch;
    if (rank_out) {
        *rank_out = coeff_rank(got->p, r);
    }
    if (r == 0) {
        if (verbose) {
            printf("[目的] batch=%d 没有收到包\n", batch);
        }
        return;
    }
    if (verbose) {
        printf("[目的] batch=%d 收到 %d 个，不再乘随机系数\n", batch, r);
    }
    plan = make_plan(code_seed, (uint32_t)batch, K);
    eq->n_recv = r;
    eq->Y = xmalloc((size_t)r * (size_t)T);
    for (j = 0; j < r; j++) {
        memcpy(eq->Y + (size_t)j * (size_t)T, got->p[j].payload, (size_t)T);
        if (trace) {
            print_coeff_payload("    [目的] 收到", got->p[j].batch_id, j,
                                got->p[j].coeff, got->p[j].payload, T, "");
        }
    }
    for (a = 0; a < plan.d; a++) {
        uint8_t *row = xmalloc((size_t)r);
        for (j = 0; j < r; j++) {
            uint8_t acc = 0;
            int m;
            for (m = 0; m < M; m++) {
                acc ^= gf_mul(plan.G[(size_t)a * M + (size_t)m], got->p[j].coeff[m]);
            }
            row[j] = acc;
        }
        batch_push_row(eq, plan.sel[a], row);
        free(row);
    }
    plan_free(&plan);
    vec_free(got);
}
