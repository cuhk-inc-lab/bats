/* 按 Yang & Yeung, Batched Sparse Codes 的 (P1) 优化外码度数分布。
   输入是一次测试里目的端转移矩阵的秩个数（秩 0..M）。
   在 x_i = (1-η)*i/N 上把 (P1) 放松成线性规划，N 默认 100。
   最大度数默认取论文里足够用的 ceil(M/η)-1。 */
#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

#include "psi_builtin.h"

#define MAX_M 64
#define MAX_SAMPLES 800
#define MAX_DEGREE 4096

/* 不给参数时用内置的那次目的端测量。写死的 Ψ 在 psi_builtin.h。 */
static const int kDefaultCounts[] = {0, 0, 0, 0, 1, 2, 11, 100, 38};
static const int kOldW[PSI_BUILTIN_COUNT] = PSI_BUILTIN_W_INIT;

static void usage(void) {
    fprintf(stderr,
            "用法: psi_opt [选项] h0 h1 ... hM\n"
            "  h_r 是秩为 r 的 batch 个数。个数一共 M+1 个，M 由此确定。\n"
            "  不给 h 时，用内置的那次目的端测量：0 0 0 0 1 2 11 100 38\n"
            "选项:\n"
            "  --eta-bar F     BP 要稳住的比例，默认 0.99。η = 1-F\n"
            "  --samples N     论文里的采样点数，默认 100\n"
            "  --q Q           有限域大小，默认 256\n"
            "  --max-degree D  度数上限。默认用论文的 ceil(M/η)-1\n"
            "  --scale S       输出整数权重的分母，默认 1000000\n"
            "  --weights-out 文件   写出度数权重，给 bats_sim --psi 用\n");
}

static double zeta_m_r(int m, int r, double q) {
    int j;
    double z = 1.0;
    if (r <= 0) {
        return 1.0;
    }
    for (j = 0; j < r; j++) {
        int e = m - j;
        if (e <= 0) {
            return 0.0;
        }
        z *= 1.0 - pow(q, (double)(-e));
    }
    return z;
}

/* I_{a,b}(x) = sum_{j=a}^{a+b-1} C(a+b-1, j) x^j (1-x)^{a+b-1-j} */
static double reg_ibeta(int a, int b, double x) {
    int n;
    int t;
    int j;
    double log_term;
    double maxlog;
    double sum;
    double logs[MAX_M + 1];
    if (b <= 0 || b > MAX_M) {
        return 0.0;
    }
    if (a <= 0) {
        return 1.0;
    }
    n = a + b - 1;
    if (x <= 0.0) {
        return 0.0;
    }
    if (x >= 1.0) {
        return 1.0;
    }
    log_term = (double)n * log(x);
    logs[0] = log_term;
    maxlog = log_term;
    for (t = 1, j = n; t < b; t++, j--) {
        log_term += log((double)j) - log((double)(n - j + 1)) + log(1.0 - x) - log(x);
        logs[t] = log_term;
        if (log_term > maxlog) {
            maxlog = log_term;
        }
    }
    if (maxlog < -700.0) {
        return 0.0;
    }
    sum = 0.0;
    for (t = 0; t < b; t++) {
        sum += exp(logs[t] - maxlog);
    }
    return exp(maxlog) * sum;
}

static int ibeta_self_test(void) {
    double a = reg_ibeta(1, 1, 0.3);
    double b = reg_ibeta(2, 1, 0.5);
    double c = reg_ibeta(1, 2, 0.25);
    double c_want = 2.0 * 0.25 * 0.75 + 0.25 * 0.25;
    if (fabs(a - 0.3) > 1e-12 || fabs(b - 0.25) > 1e-12 || fabs(c - c_want) > 1e-12) {
        fprintf(stderr, "不完全贝塔自检失败\n");
        return 0;
    }
    return 1;
}

static void fill_hbar(const double *h, int M, double q, double *hbar) {
    int r;
    int i;
    for (r = 1; r <= M; r++) {
        double s = 0.0;
        for (i = r; i <= M; i++) {
            double z = zeta_m_r(i, r, q);
            s += z / pow(q, (double)(i - r)) * h[i];
        }
        hbar[r] = s;
    }
}

static double hbar_tail(const double *hbar, int M, int r0) {
    int r;
    double s = 0.0;
    if (r0 < 1) {
        r0 = 1;
    }
    for (r = r0; r <= M; r++) {
        s += hbar[r];
    }
    return s;
}

static double hprime_direct(const double *h, int M, double q, int r) {
    int i;
    double s = 0.0;
    for (i = r; i <= M; i++) {
        s += zeta_m_r(i, r, q) * h[i];
    }
    return s;
}

/* (30)：c_d(x) 使 Ω = sum_d Ψ_d c_d(x) */
static double cost_of(int d, double x, const double *hbar, int M) {
    int r;
    int rmax;
    double acc = 0.0;
    rmax = d - 1;
    if (rmax > M) {
        rmax = M;
    }
    for (r = 1; r <= rmax; r++) {
        acc += hbar[r] * reg_ibeta(d - r, r, x);
    }
    if (d <= M) {
        acc += hbar_tail(hbar, M, d);
    }
    return (double)d * acc;
}

static double omega_of(const double *psi, int D, double x, const double *hbar, int M) {
    int d;
    double s = 0.0;
    for (d = 1; d <= D; d++) {
        if (psi[d] == 0.0) {
            continue;
        }
        s += psi[d] * cost_of(d, x, hbar, M);
    }
    return s;
}

static double theta_on_grid(const double *psi, int D, const double *xs, int n, const double *hbar,
                            int M) {
    int i;
    double best = 1e300;
    for (i = 0; i < n; i++) {
        double om = omega_of(psi, D, xs[i], hbar, M);
        double gap = -log(1.0 - xs[i]);
        double t;
        if (gap <= 0.0) {
            continue;
        }
        t = om / gap;
        if (t < best) {
            best = t;
        }
    }
    if (best > 1e299) {
        return 0.0;
    }
    return best;
}

static int invert_basis(const double *A, int m, int n, const int *basic, const double *rhs,
                        double *x, double *Binv) {
    int stride = 2 * m + 1;
    double *aug = (double *)malloc((size_t)m * (size_t)stride * sizeof(double));
    int i;
    int j;
    int k;
    if (!aug) {
        return 0;
    }
    for (i = 0; i < m; i++) {
        for (j = 0; j < m; j++) {
            aug[i * stride + j] = A[i * n + basic[j]];
        }
        for (j = 0; j < m; j++) {
            aug[i * stride + m + j] = (i == j) ? 1.0 : 0.0;
        }
        aug[i * stride + 2 * m] = rhs[i];
    }
    for (k = 0; k < m; k++) {
        int piv = k;
        double best = fabs(aug[k * stride + k]);
        double diag;
        for (i = k + 1; i < m; i++) {
            double v = fabs(aug[i * stride + k]);
            if (v > best) {
                best = v;
                piv = i;
            }
        }
        if (best < 1e-14) {
            free(aug);
            return 0;
        }
        if (piv != k) {
            for (j = k; j < stride; j++) {
                double tmp = aug[k * stride + j];
                aug[k * stride + j] = aug[piv * stride + j];
                aug[piv * stride + j] = tmp;
            }
        }
        diag = aug[k * stride + k];
        for (j = k; j < stride; j++) {
            aug[k * stride + j] /= diag;
        }
        for (i = 0; i < m; i++) {
            double f;
            if (i == k) {
                continue;
            }
            f = aug[i * stride + k];
            if (f == 0.0) {
                continue;
            }
            for (j = k; j < stride; j++) {
                aug[i * stride + j] -= f * aug[k * stride + j];
            }
        }
    }
    for (i = 0; i < m; i++) {
        x[i] = aug[i * stride + 2 * m];
        for (j = 0; j < m; j++) {
            Binv[i * m + j] = aug[i * stride + m + j];
        }
    }
    free(aug);
    return 1;
}

static int simplex(const double *A, const double *rhs, const double *cost, int m, int n,
                   const int *basic0, double *psi, int D, double *theta, int *pivots_out) {
    int *basic = (int *)malloc((size_t)m * sizeof(int));
    int *is_basic = (int *)calloc((size_t)n, sizeof(int));
    double *x = (double *)malloc((size_t)m * sizeof(double));
    double *Binv = (double *)malloc((size_t)m * (size_t)m * sizeof(double));
    double *pi = (double *)malloc((size_t)m * sizeof(double));
    double *abar = (double *)malloc((size_t)m * sizeof(double));
    double best_theta = -1.0;
    double *best_psi = (double *)calloc((size_t)(D + 1), sizeof(double));
    int pivots = 0;
    int ok = 0;
    int stall = 0;
    int bland = 0;
    int i;
    if (!basic || !is_basic || !x || !Binv || !pi || !abar || !best_psi) {
        fprintf(stderr, "内存不足\n");
        goto done;
    }
    memcpy(basic, basic0, (size_t)m * sizeof(int));
    for (i = 0; i < m; i++) {
        is_basic[basic[i]] = 1;
    }
    while (pivots < 40000) {
        double obj;
        int enter = -1;
        double best_cbar = 0.0;
        int leave = -1;
        double best_ratio = 0.0;
        int leave_var = 0;
        double resid = 0.0;
        int feasible = 1;
        if (!invert_basis(A, m, n, basic, rhs, x, Binv)) {
            fprintf(stderr, "基矩阵奇异，停在第 %d 次换基\n", pivots);
            break;
        }
        obj = 0.0;
        for (i = 0; i < m; i++) {
            double s = 0.0;
            int k;
            for (k = 0; k < m; k++) {
                s += A[i * n + basic[k]] * x[k];
            }
            resid += fabs(s - rhs[i]);
            obj += cost[basic[i]] * x[i];
            if (x[i] < -1e-7) {
                feasible = 0;
            }
        }
        if (resid > 1e-5) {
            fprintf(stderr, "基的残差 %.3g 过大\n", resid);
            break;
        }
        if (feasible && obj > best_theta + 1e-12) {
            memset(best_psi, 0, (size_t)(D + 1) * sizeof(double));
            best_theta = obj;
            for (i = 0; i < m; i++) {
                int v = basic[i];
                if (v < D && x[i] > 0.0) {
                    best_psi[v + 1] = x[i];
                }
            }
            stall = 0;
            bland = 0;
        } else if (feasible) {
            stall++;
        }
        for (i = 0; i < m; i++) {
            int k;
            double pik = 0.0;
            for (k = 0; k < m; k++) {
                pik += cost[basic[k]] * Binv[k * m + i];
            }
            pi[i] = pik;
        }
        best_cbar = bland ? 0.0 : 1e-9;
        enter = -1;
        for (i = 0; i < n; i++) {
            int r;
            double dot = 0.0;
            double cbar;
            if (is_basic[i]) {
                continue;
            }
            for (r = 0; r < m; r++) {
                dot += pi[r] * A[r * n + i];
            }
            cbar = cost[i] - dot;
            if (bland) {
                if (cbar > 1e-9) {
                    enter = i;
                    best_cbar = cbar;
                    break;
                }
            } else if (cbar > best_cbar) {
                best_cbar = cbar;
                enter = i;
            }
        }
        if (enter < 0) {
            ok = 1;
            break;
        }
        for (i = 0; i < m; i++) {
            int r;
            double s = 0.0;
            for (r = 0; r < m; r++) {
                s += Binv[i * m + r] * A[r * n + enter];
            }
            abar[i] = s;
        }
        leave = -1;
        best_ratio = 0.0;
        leave_var = n + 1;
        for (i = 0; i < m; i++) {
            double ratio;
            int bvar;
            if (abar[i] <= 1e-12) {
                continue;
            }
            ratio = x[i] / abar[i];
            if (ratio < -1e-8) {
                continue;
            }
            if (ratio < 0.0) {
                ratio = 0.0;
            }
            bvar = basic[i];
            if (leave < 0 || ratio < best_ratio - 1e-14 ||
                (fabs(ratio - best_ratio) <= 1e-14 && bvar < leave_var)) {
                best_ratio = ratio;
                leave = i;
                leave_var = bvar;
            }
        }
        if (leave < 0) {
            fprintf(stderr, "线性规划无界\n");
            break;
        }
        is_basic[basic[leave]] = 0;
        basic[leave] = enter;
        is_basic[enter] = 1;
        pivots++;
        if (stall > 80) {
            bland = 1;
        }
        if (pivots % 500 == 0) {
            fprintf(stderr, "单纯形已换基 %d 次，当前 θ=%.6f\n", pivots, best_theta);
        }
    }
    if (pivots >= 40000) {
        fprintf(stderr, "达到换基上限，采用目前最好的可行解\n");
    }
    if (best_theta >= 0.0) {
        double sum = 0.0;
        int d;
        for (d = 1; d <= D; d++) {
            if (best_psi[d] < 0.0) {
                best_psi[d] = 0.0;
            }
            sum += best_psi[d];
        }
        if (sum > 0.0) {
            for (d = 1; d <= D; d++) {
                psi[d] = best_psi[d] / sum;
            }
            *theta = best_theta;
            ok = 1;
        }
    }
done:
    if (pivots_out) {
        *pivots_out = pivots;
    }
    free(basic);
    free(is_basic);
    free(x);
    free(Binv);
    free(pi);
    free(abar);
    free(best_psi);
    return ok;
}

static int solve_p1(const double *xs, int N, int D, int M, const double *hbar, double *psi,
                    double *theta, int *pivots) {
    int m = N + 1;
    int n = D + 1 + N;
    int theta_col = D;
    double *A = (double *)calloc((size_t)m * (size_t)n, sizeof(double));
    double *rhs = (double *)calloc((size_t)m, sizeof(double));
    double *cost = (double *)calloc((size_t)n, sizeof(double));
    int *basic = (int *)malloc((size_t)m * sizeof(int));
    double *cd = (double *)malloc((size_t)N * (size_t)(D + 1) * sizeof(double));
    int i;
    int d;
    int rc;
    if (!A || !rhs || !cost || !basic || !cd) {
        fprintf(stderr, "内存不足\n");
        free(A);
        free(rhs);
        free(cost);
        free(basic);
        free(cd);
        return 0;
    }
    for (i = 0; i < N; i++) {
        for (d = 1; d <= D; d++) {
            cd[(size_t)i * (size_t)(D + 1) + (size_t)d] = cost_of(d, xs[i], hbar, M);
        }
    }
    rhs[0] = 1.0;
    for (d = 1; d <= D; d++) {
        A[d - 1] = 1.0;
    }
    basic[0] = 0;
    for (i = 0; i < N; i++) {
        int row = i + 1;
        double gap = -log(1.0 - xs[i]);
        double scale = gap;
        int col;
        for (d = 1; d <= D; d++) {
            double c = cd[(size_t)i * (size_t)(D + 1) + (size_t)d];
            if (c > scale) {
                scale = c;
            }
        }
        if (scale < 1.0) {
            scale = 1.0;
        }
        for (d = 1; d <= D; d++) {
            A[row * n + (d - 1)] = cd[(size_t)i * (size_t)(D + 1) + (size_t)d] / scale;
        }
        A[row * n + theta_col] = -gap / scale;
        col = D + 1 + i;
        A[row * n + col] = -1.0;
        rhs[row] = 0.0;
        basic[row] = col;
    }
    cost[theta_col] = 1.0;
    rc = simplex(A, rhs, cost, m, n, basic, psi, D, theta, pivots);
    free(cd);
    free(A);
    free(rhs);
    free(cost);
    free(basic);
    return rc;
}

static int write_weights(const char *path, const int *w, int D) {
    FILE *f = fopen(path, "wb");
    int d;
    if (!f) {
        fprintf(stderr, "写权重失败 %s\n", path);
        return 0;
    }
    fprintf(f, "# BATS 度数权重，度数 1..D。给 bats_sim --psi 用。\n");
    fprintf(f, "%d\n", D);
    for (d = 1; d <= D; d++) {
        fprintf(f, "%d%c", w[d], d == D ? '\n' : ' ');
    }
    fclose(f);
    return 1;
}

static int read_counts_file(const char *path, int *counts, int cap, int *n_out) {
    FILE *f = fopen(path, "rb");
    int n = 0;
    int v;
    int c;
    if (!f) {
        fprintf(stderr, "打不开 %s\n", path);
        return 0;
    }
    while (n < cap) {
        c = fgetc(f);
        if (c == EOF) {
            break;
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
        if (fscanf(f, "%d", &v) != 1) {
            fprintf(stderr, "%s 里有非整数\n", path);
            fclose(f);
            return 0;
        }
        if (v < 0) {
            fprintf(stderr, "秩个数不能是负数\n");
            fclose(f);
            return 0;
        }
        counts[n++] = v;
    }
    fclose(f);
    if (n < 2) {
        fprintf(stderr, "%s 里的秩个数不够\n", path);
        return 0;
    }
    *n_out = n;
    return 1;
}

int psi_optimize(const int *counts, int ncounts, int max_degree, double eta_bar, int samples,
                 double q, int scale, int *weights_out, int weights_cap, int *D_out,
                 const char *weights_path) {
    int M;
    int D;
    int Dpaper;
    int N = samples;
    double eta;
    double h[MAX_M + 1];
    double hbar[MAX_M + 1];
    double *psi = NULL;
    double *xs = NULL;
    double *fine = NULL;
    int *weights = NULL;
    int i;
    int d;
    int round;
    int pivots = 0;
    long long total = 0;
    double hmean = 0.0;
    double theta = 0.0;
    double theta_fine;
    double theta_int;
    double sum_psi;
    double mean_deg;
    double omega0;
    double lemma_err = 0.0;
    int rstar = 0;
    int rc = 0;
    if (!ibeta_self_test()) {
        return 0;
    }
    if (!counts || ncounts < 2 || ncounts > MAX_M + 1) {
        fprintf(stderr, "需要 2 到 %d 个秩个数（秩 0 到 M）\n", MAX_M + 1);
        return 0;
    }
    if (eta_bar <= 0.0 || eta_bar >= 1.0 || N < 2 || N > MAX_SAMPLES || q <= 2.0 || scale < 1) {
        fprintf(stderr, "优化参数不合法\n");
        return 0;
    }
    M = ncounts - 1;
    eta = 1.0 - eta_bar;
    Dpaper = (int)ceil((double)M / eta) - 1;
    if (Dpaper < 1) {
        Dpaper = 1;
    }
    if (Dpaper > MAX_DEGREE) {
        Dpaper = MAX_DEGREE;
    }
    D = Dpaper;
    if (max_degree > 0 && max_degree < D) {
        D = max_degree;
    }
    if (weights_out && weights_cap < D + 1) {
        fprintf(stderr, "权重缓冲区太小\n");
        return 0;
    }
    for (i = 0; i <= M; i++) {
        if (counts[i] < 0) {
            fprintf(stderr, "秩个数不能是负数\n");
            return 0;
        }
        total += counts[i];
    }
    if (total <= 0) {
        fprintf(stderr, "秩个数全是 0\n");
        return 0;
    }
    for (i = 0; i <= M; i++) {
        h[i] = (double)counts[i] / (double)total;
        hmean += (double)i * h[i];
        if (counts[i] > 0) {
            rstar = i;
        }
    }
    fill_hbar(h, M, q, hbar);
    for (i = 1; i <= M; i++) {
        double diff = fabs(hprime_direct(h, M, q, i) - hbar_tail(hbar, M, i));
        if (diff > lemma_err) {
            lemma_err = diff;
        }
    }
    if (lemma_err > 1e-6) {
        fprintf(stderr, "有效秩恒等式对不上，最大偏差 %.3g\n", lemma_err);
        return 0;
    }
    psi = (double *)calloc((size_t)(D + 1), sizeof(double));
    xs = (double *)malloc((size_t)MAX_SAMPLES * sizeof(double));
    fine = (double *)malloc(4000 * sizeof(double));
    weights = (int *)calloc((size_t)(D + 1), sizeof(int));
    if (!psi || !xs || !fine || !weights) {
        fprintf(stderr, "内存不足\n");
        goto done;
    }
    printf("q=%.0f  M=%d  样本 batch=%lld  平均秩 h_bar=%.4f\n", q, M, total, hmean);
    printf("有效秩 ħ_r（随机 G 乘上转移矩阵之后，BP 用的秩）:");
    for (i = 1; i <= M; i++) {
        printf(" %.4f", hbar[i]);
    }
    printf("\n");
    printf("η_bar=%.4f  论文足够的最大度数 D=%d", eta_bar, Dpaper);
    if (D != Dpaper) {
        printf("，这次限制到 %d", D);
    }
    printf("  采样 N=%d\n", N);
    for (i = 0; i < N; i++) {
        xs[i] = eta_bar * (double)(i + 1) / (double)N;
    }
    for (round = 0; round < 6; round++) {
        int nfine = 4000;
        int added = 0;
        double worst = 1e300;
        printf("第 %d 轮线性规划，约束点 %d 个\n", round + 1, N);
        fflush(stdout);
        if (!solve_p1(xs, N, D, M, hbar, psi, &theta, &pivots)) {
            fprintf(stderr, "线性规划没有可行解\n");
            goto done;
        }
        printf("换基 %d 次，采样点上的 θ=%.6f\n", pivots, theta);
        for (i = 0; i < nfine; i++) {
            fine[i] = eta_bar * (double)(i + 1) / (double)nfine;
        }
        theta_fine = theta_on_grid(psi, D, fine, nfine, hbar, M);
        printf("细网格上真正能保证的 θ=%.6f\n", theta_fine);
        if (theta_fine >= theta - 1e-4 * (1.0 + fabs(theta))) {
            theta = theta_fine;
            break;
        }
        if (N + 8 > MAX_SAMPLES) {
            theta = theta_fine;
            fprintf(stderr, "采样点已满，采用细网格上的 θ\n");
            break;
        }
        for (added = 0; added < 8 && N < MAX_SAMPLES; added++) {
            int worst_i = -1;
            worst = 1e300;
            for (i = 0; i < nfine; i++) {
                double om = omega_of(psi, D, fine[i], hbar, M);
                double t = om / -log(1.0 - fine[i]);
                int s;
                int dup = 0;
                if (t >= worst) {
                    continue;
                }
                for (s = 0; s < N; s++) {
                    if (fabs(xs[s] - fine[i]) < 1e-12) {
                        dup = 1;
                        break;
                    }
                }
                if (!dup) {
                    worst = t;
                    worst_i = i;
                }
            }
            if (worst_i < 0 || worst >= theta - 1e-4 * (1.0 + fabs(theta))) {
                break;
            }
            xs[N++] = fine[worst_i];
        }
        if (added == 0) {
            theta = theta_fine;
            break;
        }
        printf("有采样点之间的缺口，补上 %d 个再解\n", added);
    }
    omega0 = 0.0;
    for (d = 1; d <= D && d <= M; d++) {
        omega0 += (double)d * psi[d] * hbar_tail(hbar, M, d);
    }
    if (omega0 < 1e-8 && rstar > 0) {
        double eps = 0.01;
        int nlow = rstar;
        printf("Ω(0)=0，按论文给度数 1..%d 补上一小块质量，否则 BP 起不来\n", rstar);
        for (d = 1; d <= D; d++) {
            psi[d] *= 1.0 - eps;
        }
        for (d = 1; d <= nlow; d++) {
            psi[d] += eps / (double)nlow;
        }
        for (i = 0; i < 4000; i++) {
            fine[i] = eta_bar * (double)(i + 1) / 4000.0;
        }
        theta = theta_on_grid(psi, D, fine, 4000, hbar, M);
    }
    sum_psi = 0.0;
    mean_deg = 0.0;
    for (d = 1; d <= D; d++) {
        sum_psi += psi[d];
        mean_deg += (double)d * psi[d];
    }
    printf("\nP1 最优 θ = %.6f\n", theta);
    printf("设计速率 η_bar*θ = %.6f 个源包/批（BP 渐近能稳住的部分）\n", eta_bar * theta);
    printf("信道容量上界 sum r*h_r = %.6f 个源包/批\n", hmean);
    printf("渐近上，BP 解出 %.2f%% 的源包大约要 K/θ = K/%.4f 个 batch\n", eta_bar * 100.0, theta);
    if (eta_bar * theta > hmean + 1e-3) {
        fprintf(stderr, "设计速率超过了容量上界，结果不可用\n");
        goto done;
    }
    printf("平均度数 = %.3f\n", mean_deg);
    printf("度数分布（只列出质量 > 1e-4 的度数）:\n");
    for (d = 1; d <= D; d++) {
        if (psi[d] > 1e-4) {
            printf("  度数 %d    %.6f\n", d, psi[d]);
        }
    }
    {
        int sumw = 0;
        int biggest = 1;
        long long acc = 0;
        for (d = 1; d <= D; d++) {
            double raw = psi[d] * (double)scale;
            int w = (int)floor(raw + 0.5);
            if (w < 0) {
                w = 0;
            }
            weights[d] = w;
            acc += w;
            if (w > weights[biggest]) {
                biggest = d;
            }
        }
        if (acc > 2000000000LL) {
            fprintf(stderr, "权重和太大\n");
            goto done;
        }
        sumw = (int)acc;
        weights[biggest] += scale - sumw;
        if (weights[biggest] < 0) {
            fprintf(stderr, "整数权重四舍五入失败\n");
            goto done;
        }
        {
            double *qpsi = (double *)calloc((size_t)(D + 1), sizeof(double));
            if (!qpsi) {
                fprintf(stderr, "内存不足\n");
                goto done;
            }
            for (d = 1; d <= D; d++) {
                qpsi[d] = (double)weights[d] / (double)scale;
            }
            for (i = 0; i < 4000; i++) {
                fine[i] = eta_bar * (double)(i + 1) / 4000.0;
            }
            theta_int = theta_on_grid(qpsi, D, fine, 4000, hbar, M);
            printf("整数权重（分母 %d）在细网格上的 θ = %.6f\n", scale, theta_int);
            printf("整数权重:");
            for (d = 1; d <= D; d++) {
                if (weights[d] > 0) {
                    printf(" %d:%d", d, weights[d]);
                }
            }
            printf("\n");
            free(qpsi);
        }
    }
    if (M == 8) {
        double *old = (double *)calloc((size_t)(D + 1), sizeof(double));
        double old_theta;
        double old_mean = 0.0;
        if (!old) {
            fprintf(stderr, "内存不足\n");
            goto done;
        }
        for (d = 1; d <= PSI_BUILTIN_COUNT && d <= D; d++) {
            old[d] = (double)kOldW[d - 1] / (double)PSI_BUILTIN_SUM;
            old_mean += (double)d * old[d];
        }
        for (i = 0; i < 4000; i++) {
            fine[i] = eta_bar * (double)(i + 1) / 4000.0;
        }
        old_theta = theta_on_grid(old, D, fine, 4000, hbar, M);
        printf("原来写死的 Ψ（平均度数 %.3f）在同一条细网格上的 θ = %.6f\n", old_mean, old_theta);
        free(old);
    }
    if (weights_path) {
        if (!write_weights(weights_path, weights, D)) {
            goto done;
        }
        printf("权重已写入 %s\n", weights_path);
    }
    if (weights_out) {
        memset(weights_out, 0, (size_t)weights_cap * sizeof(int));
        memcpy(weights_out, weights, (size_t)(D + 1) * sizeof(int));
    }
    if (D_out) {
        *D_out = D;
    }
    printf("和为 1 的检查: %.6f\n", sum_psi);
    rc = 1;
done:
    free(psi);
    free(xs);
    free(fine);
    free(weights);
    return rc;
}

#ifdef PSI_OPT_MAIN
int main(int argc, char **argv) {
    int counts[MAX_M + 1];
    int ncounts = 0;
    int N = 100;
    int scale = 1000000;
    int max_degree = 0;
    double eta_bar = 0.99;
    double q = 256.0;
    const char *out_path = NULL;
    int arg;
#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif
    for (arg = 1; arg < argc; arg++) {
        if (strcmp(argv[arg], "--eta-bar") == 0 && arg + 1 < argc) {
            eta_bar = atof(argv[++arg]);
        } else if (strcmp(argv[arg], "--samples") == 0 && arg + 1 < argc) {
            N = atoi(argv[++arg]);
        } else if (strcmp(argv[arg], "--q") == 0 && arg + 1 < argc) {
            q = atof(argv[++arg]);
        } else if (strcmp(argv[arg], "--max-degree") == 0 && arg + 1 < argc) {
            max_degree = atoi(argv[++arg]);
        } else if (strcmp(argv[arg], "--scale") == 0 && arg + 1 < argc) {
            scale = atoi(argv[++arg]);
        } else if (strcmp(argv[arg], "--weights-out") == 0 && arg + 1 < argc) {
            out_path = argv[++arg];
        } else if (strcmp(argv[arg], "--help") == 0) {
            usage();
            return 0;
        } else if (argv[arg][0] == '-') {
            usage();
            return 1;
        } else {
            char *end = NULL;
            long v = strtol(argv[arg], &end, 10);
            if (end != argv[arg] && *end == '\0') {
                if (v < 0 || ncounts > MAX_M) {
                    fprintf(stderr, "秩个数不合法\n");
                    return 1;
                }
                counts[ncounts++] = (int)v;
            } else {
                if (!read_counts_file(argv[arg], counts, MAX_M + 1, &ncounts)) {
                    return 1;
                }
            }
        }
    }
    if (ncounts == 0) {
        ncounts = (int)(sizeof(kDefaultCounts) / sizeof(kDefaultCounts[0]));
        memcpy(counts, kDefaultCounts, sizeof(kDefaultCounts));
        printf("未给秩个数，用内置测量（目的端，刚好能还原的前 152 个 batch）\n");
    }
    return psi_optimize(counts, ncounts, max_degree, eta_bar, N, q, scale, NULL, 0, NULL,
                        out_path)
               ? 0
               : 1;
}
#endif
