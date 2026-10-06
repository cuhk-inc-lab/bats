/* 外码：进程里一份 Ψ，按它抽度数并生成 G。G 不进包。发出包在 sim.c。 */
#include "bats_internal.h"
#include "psi_builtin.h"

typedef char psi_builtin_matches_m[(PSI_BUILTIN_COUNT == M) ? 1 : -1];

const int PSI_DEG[PSI_BUILTIN_COUNT] = {1, 2, 3, 4, 5, 6, 7, 8};
const int PSI_W[PSI_BUILTIN_COUNT] = PSI_BUILTIN_W_INIT;
const int PSI_SUM = PSI_BUILTIN_SUM;

#define PSI_CAP 4096

static int psi_w[PSI_CAP + 1];
static int psi_deg_max;
static int psi_sum;
static int psi_ready;
static int psi_custom;
static int psi_from_opt;

static int psi_next_int(FILE *f, int *out) {
    int c;
    for (;;) {
        c = fgetc(f);
        if (c == EOF) {
            return 0;
        }
        if (c == '#') {
            while ((c = fgetc(f)) != EOF && c != '\n') {
            }
            continue;
        }
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            continue;
        }
        ungetc(c, f);
        return fscanf(f, "%d", out) == 1;
    }
}

static void psi_use_builtin(void) {
    int i;
    memset(psi_w, 0, sizeof(psi_w));
    for (i = 0; i < M; i++) {
        psi_w[PSI_DEG[i]] = PSI_W[i];
    }
    psi_deg_max = M;
    psi_sum = PSI_SUM;
    psi_ready = 1;
}

static void psi_ensure(void) {
    if (!psi_ready) {
        psi_use_builtin();
    }
}

static void psi_store(const int *weights, int max_degree, int sum, int from_opt) {
    memset(psi_w, 0, sizeof(psi_w));
    memcpy(psi_w, weights, ((size_t)max_degree + 1) * sizeof(int));
    psi_deg_max = max_degree;
    psi_sum = sum;
    psi_ready = 1;
    psi_custom = 1;
    psi_from_opt = from_opt;
}

int psi_is_custom(void) {
    return psi_custom;
}

int psi_load(const char *path) {
    FILE *f = fopen(path, "rb");
    int D;
    int d;
    int *w;
    int sum = 0;
    if (!f) {
        fprintf(stderr, "打不开度数分布 %s\n", path);
        return 0;
    }
    if (!psi_next_int(f, &D) || D < 1 || D > PSI_CAP) {
        fprintf(stderr, "%s 的最大度数不合法\n", path);
        fclose(f);
        return 0;
    }
    w = (int *)calloc((size_t)D + 1, sizeof(int));
    if (!w) {
        fprintf(stderr, "内存不足\n");
        fclose(f);
        return 0;
    }
    for (d = 1; d <= D; d++) {
        int weight;
        if (!psi_next_int(f, &weight) || weight < 0) {
            fprintf(stderr, "%s 的权重不合法\n", path);
            free(w);
            fclose(f);
            return 0;
        }
        w[d] = weight;
        if (sum > 2000000000 - weight) {
            fprintf(stderr, "%s 的权重和太大\n", path);
            free(w);
            fclose(f);
            return 0;
        }
        sum += weight;
    }
    fclose(f);
    if (sum <= 0) {
        fprintf(stderr, "%s 的权重全是 0\n", path);
        free(w);
        return 0;
    }
    psi_store(w, D, sum, 0);
    free(w);
    return 1;
}

int psi_set(const int *weights, int max_degree) {
    int d;
    int sum = 0;
    if (!weights || max_degree < 1 || max_degree > PSI_CAP) {
        fprintf(stderr, "度数分布不合法\n");
        return 0;
    }
    for (d = 1; d <= max_degree; d++) {
        if (weights[d] < 0) {
            fprintf(stderr, "度数分布不合法\n");
            return 0;
        }
        if (sum > 2000000000 - weights[d]) {
            fprintf(stderr, "权重和太大\n");
            return 0;
        }
        sum += weights[d];
    }
    if (sum <= 0) {
        fprintf(stderr, "度数分布的权重全是 0\n");
        return 0;
    }
    psi_store(weights, max_degree, sum, 1);
    return 1;
}

void psi_print(FILE *fp) {
    int d;
    psi_ensure();
    fprintf(fp, "%s（度数:权重，分母 %d）:",
            !psi_custom ? "写死的 Ψ" : (psi_from_opt ? "优化得到的 Ψ" : "载入的 Ψ"), psi_sum);
    for (d = 1; d <= psi_deg_max; d++) {
        if (psi_w[d] > 0) {
            fprintf(fp, " %d:%d", d, psi_w[d]);
        }
    }
    fprintf(fp, "\n");
}

static int sample_degree(Rng *rng) {
    int x;
    int acc = 0;
    int d;
    psi_ensure();
    x = (int)rng_below(rng, (uint32_t)psi_sum);
    for (d = 1; d <= psi_deg_max; d++) {
        acc += psi_w[d];
        if (x < acc) {
            return d;
        }
    }
    return psi_deg_max;
}

/* 线上 Ψ 的度数会超过 32。缓冲放到线程本地，避免每个 batch 都 malloc。 */
enum { PLAN_D = 512, PLAN_K = 512 };
static _Thread_local int plan_tls_busy;
static _Thread_local int plan_tls_sel[PLAN_D];
static _Thread_local uint8_t plan_tls_G[PLAN_D * M];
static _Thread_local int plan_tls_pool[PLAN_K];

Plan make_plan(uint64_t seed, uint32_t batch_id, int K) {
    Plan plan;
    Rng rng = rng_for_batch(seed, batch_id);
    int d = sample_degree(&rng);
    int *pool;
    int i;
    memset(&plan, 0, sizeof(plan));
    if (d > K) {
        d = K;
    }
    plan.d = d;
    if (d <= PLAN_D && K <= PLAN_K && !plan_tls_busy) {
        plan.sel = plan_tls_sel;
        plan.G = plan_tls_G;
        pool = plan_tls_pool;
        plan.heap = 0;
        plan_tls_busy = 1;
    } else {
        plan.sel = xmalloc((size_t)d * sizeof(int));
        plan.G = xmalloc((size_t)d * (size_t)M);
        pool = xmalloc((size_t)K * sizeof(int));
        plan.heap = 1;
    }
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
    if (plan.heap) {
        free(pool);
    }
    return plan;
}

void plan_free(Plan *plan) {
    if (plan->heap) {
        free(plan->sel);
        free(plan->G);
    } else {
        plan_tls_busy = 0;
    }
    plan->sel = NULL;
    plan->G = NULL;
    plan->heap = 0;
}

void psi_reset(void) {
    psi_use_builtin();
    psi_custom = 0;
    psi_from_opt = 0;
}

int psi_weight_sum(void) {
    psi_ensure();
    return psi_sum;
}

long long psi_degree_moment(void) {
    int d;
    long long deg_sum = 0;
    psi_ensure();
    for (d = 1; d <= psi_deg_max; d++) {
        deg_sum += (long long)d * (long long)psi_w[d];
    }
    return deg_sum;
}
