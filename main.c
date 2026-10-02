/* 自检、小例子、开销测量，以及按文件还原。 */
#include "bats.h"

#ifdef _WIN32
#include <windows.h>
#endif

enum { K_MEASURE = 128 };
static const uint64_t kCodeSeed = 0xBA75C0DEULL;
static const uint64_t kLossSeed = 0x1055EE01ULL;
static const uint64_t kRecodeSeed = 0xBEC0DE01ULL;

static int run_silent(const uint8_t *src, int K, int T, uint64_t seed, int loss, int n_batches) {
    uint8_t *dst = xmalloc((size_t)K * (size_t)T);
    SimResult r = simulate(src, K, T, seed, seed ^ 0x1111ULL, seed ^ 0x2222ULL, loss,
                           n_batches, 0, -1, dst);
    free(dst);
    return r.ok;
}

static int run_example(void) {
    enum { K = 4, T = 1, N_BATCH = 6 };
    uint8_t src[K] = {0x42, 0x41, 0x54, 0x53}; /* B A T S */
    uint8_t dst[K];
    SimResult r;
    int i;
    int same;
    printf("===== 小例子：每个包 1 字节 =====\n");
    printf("q=256  M=%d  节点=源、中继1、中继2、目的  丢包=%d%%  无预编码\n", M, LOSS_PERCENT);
    psi_print(stdout);
    printf("batch_id 只是箱子编号。下面展开全部 %d 个 batch 的系数、payload 和丢包后的 recode。\n",
           N_BATCH);
    r = simulate(src, K, T, kCodeSeed, kLossSeed, kRecodeSeed, LOSS_PERCENT,
                 N_BATCH, 1, -1, dst);
    printf("--- 解回来的字节 ---\n");
    same = 1;
    for (i = 0; i < K; i++) {
        int printable = dst[i] >= 32 && dst[i] < 127;
        printf("  [%d] 原文=%02X 解回=%02X", i, src[i], dst[i]);
        if (printable) {
            printf(" '%c'", dst[i]);
        }
        if (dst[i] != src[i]) {
            printf(" 不一致");
            same = 0;
        }
        printf("\n");
    }
    printf("BP 解出 %d 个源包，inactivation %d 个\n", r.bp_solved, r.inactivated);
    printf("三跳 发送/丢弃/送达: %d/%d/%d , %d/%d/%d , %d/%d/%d\n",
           r.link[0].sent, r.link[0].dropped, r.link[0].delivered,
           r.link[1].sent, r.link[1].dropped, r.link[1].delivered,
           r.link[2].sent, r.link[2].dropped, r.link[2].delivered);
    if (!same || !r.ok) {
        printf("小例子未能还原\n");
        return 0;
    }
    printf("小例子与原文一致\n");
    return 1;
}

static int apply_psi(const int *hist, int measure_k, const char *weights_path) {
    enum { WCAP = 4097 };
    int weights[WCAP];
    int D = 0;
    if (!psi_optimize(hist, M + 1, measure_k, 0.99, 100, 256.0, 1000000, weights, WCAP, &D,
                      weights_path)) {
        return 0;
    }
    return psi_set(weights, D);
}

static int file_try(const uint8_t *src, uint8_t *dst, int K, int T) {
    int n_batches = batches_needed(K);
    SimResult r;
    printf("q=256 M=%d 批次数=%d 丢包=%d%% 四节点内存队列\n", M, n_batches, LOSS_PERCENT);
    r = simulate(src, K, T, kCodeSeed, kLossSeed, kRecodeSeed, LOSS_PERCENT, n_batches, 0, -1,
                 dst);
    printf("三跳 发送/丢弃/送达: %d/%d/%d , %d/%d/%d , %d/%d/%d\n",
           r.link[0].sent, r.link[0].dropped, r.link[0].delivered,
           r.link[1].sent, r.link[1].dropped, r.link[1].delivered,
           r.link[2].sent, r.link[2].dropped, r.link[2].delivered);
    printf("BP 解出 %d 个源包，inactivation %d 个\n", r.bp_solved, r.inactivated);
    return r.ok && memcmp(src, dst, (size_t)K * (size_t)T) == 0;
}

static int run_file(const char *in_path, const char *out_path, const int *hist, int optimize) {
    FILE *f;
    long len;
    uint8_t *file = NULL;
    uint8_t *src = NULL;
    uint8_t *dst = NULL;
    int T = FILE_T;
    int K;
    int matched;
    size_t nread;
    f = fopen(in_path, "rb");
    if (!f) {
        fprintf(stderr, "打不开输入文件 %s\n", in_path);
        return 0;
    }
    if (fseek(f, 0, SEEK_END) != 0) {
        fclose(f);
        fprintf(stderr, "读输入失败\n");
        return 0;
    }
    len = ftell(f);
    if (len < 0) {
        fclose(f);
        fprintf(stderr, "读输入失败\n");
        return 0;
    }
    if (fseek(f, 0, SEEK_SET) != 0) {
        fclose(f);
        fprintf(stderr, "读输入失败\n");
        return 0;
    }
    file = xmalloc((size_t)len + 1);
    nread = fread(file, 1, (size_t)len, f);
    fclose(f);
    if (nread != (size_t)len) {
        free(file);
        fprintf(stderr, "读输入失败\n");
        return 0;
    }
    if (len == 0) {
        FILE *o = fopen(out_path, "wb");
        free(file);
        if (!o) {
            fprintf(stderr, "写输出失败\n");
            return 0;
        }
        fclose(o);
        printf("空文件，视为一致\n");
        return 1;
    }
    K = (int)((len + T - 1) / T);
    if (optimize && hist && K != K_MEASURE) {
        printf("===== 文件有 %d 个源包，按这个 K 重新求 Ψ =====\n", K);
        if (!apply_psi(hist, K, NULL)) {
            free(file);
            return 0;
        }
    }
    src = xmalloc((size_t)K * (size_t)T);
    dst = xmalloc((size_t)K * (size_t)T);
    memset(src, 0, (size_t)K * (size_t)T);
    memcpy(src, file, (size_t)len);
    printf("===== 文件 =====\n");
    printf("输入 %s  长度 %ld 字节  切成 K=%d 个源包  每包 T=%d 字节（末包用 0 补齐）\n",
           in_path, len, K, T);
    matched = file_try(src, dst, K, T);
    if (!matched && optimize && psi_is_custom()) {
        fprintf(stderr, "这份 Ψ 没有还原文件，改回写死的分布再试一次\n");
        psi_reset();
        matched = file_try(src, dst, K, T);
    }
    if (!matched) {
        fprintf(stderr, "还原结果与原文不一致\n");
        free(file);
        free(src);
        free(dst);
        return 0;
    }
    {
        FILE *o = fopen(out_path, "wb");
        if (!o) {
            fprintf(stderr, "写输出失败 %s\n", out_path);
            free(file);
            free(src);
            free(dst);
            return 0;
        }
        if (fwrite(dst, 1, (size_t)len, o) != (size_t)len) {
            fprintf(stderr, "写输出失败\n");
            fclose(o);
            free(file);
            free(src);
            free(dst);
            return 0;
        }
        fclose(o);
    }
    if (memcmp(file, dst, (size_t)len) != 0) {
        fprintf(stderr, "写出的文件与原文不一致\n");
        free(file);
        free(src);
        free(dst);
        return 0;
    }
    printf("还原文件 %s 与原文一致\n", out_path);
    free(file);
    free(src);
    free(dst);
    return 1;
}

static void free_buckets(Vec *buckets, int n) {
    int i;
    if (!buckets) {
        return;
    }
    for (i = 0; i < n; i++) {
        vec_free(&buckets[i]);
    }
    free(buckets);
}

static void print_rank_row(const char *name, const int *ranks, int n) {
    int hist[M + 1];
    int i;
    long long sum = 0;
    memset(hist, 0, sizeof(hist));
    for (i = 0; i < n; i++) {
        int rk = ranks[i];
        if (rk < 0) {
            rk = 0;
        }
        if (rk > M) {
            rk = M;
        }
        hist[rk]++;
        sum += ranks[i];
    }
    printf("  %-6s", name);
    for (i = 0; i <= M; i++) {
        printf(" %4d", hist[i]);
    }
    if (n > 0) {
        printf("   %lld.%lld\n", sum / n, (sum * 10 / n) % 10);
    } else {
        printf("\n");
    }
}

static int try_prefix(const uint8_t *src, int K, int T, uint64_t code_seed, const Vec *dest,
                      int n, int *bp, int *inact) {
    BatchEq *eqs;
    uint8_t *dst;
    int ok;
    int b;
    if (n <= 0) {
        if (bp) {
            *bp = 0;
        }
        if (inact) {
            *inact = 0;
        }
        return 0;
    }
    eqs = xmalloc((size_t)n * sizeof(BatchEq));
    dst = xmalloc((size_t)K * (size_t)T);
    memset(eqs, 0, (size_t)n * sizeof(BatchEq));
    for (b = 0; b < n; b++) {
        eqs_from_packets(&eqs[b], &dest[b], code_seed, b, K, T);
    }
    ok = decode(eqs, n, K, T, dst, bp, inact);
    if (ok) {
        ok = memcmp(src, dst, (size_t)K * (size_t)T) == 0;
    }
    for (b = 0; b < n; b++) {
        batch_free(&eqs[b]);
    }
    free(eqs);
    free(dst);
    return ok;
}

/* K=128、T=32、同一条 5% 线网上，找最少 batch，并给出秩分布和开销。
   hist_out 非空时写入全部发出 batch 在目的端的秩个数，长度 M+1。 */
static int run_measure(int *hist_out) {
    enum { K = 128, T = 32 };
    uint8_t *src;
    int n_cap;
    int *rank_r1;
    int *rank_r2;
    int *rank_dst;
    Vec *dest = NULL;
    int lo;
    int hi;
    int n_min;
    int bp = 0;
    int inact = 0;
    int i;
    int ok = 0;
    long long rank_sum = 0;
    const uint64_t code_seed = kCodeSeed;
    src = xmalloc((size_t)K * (size_t)T);
    for (i = 0; i < K * T; i++) {
        src[i] = (uint8_t)((i * 17 + 3) & 255);
    }
    n_cap = batches_needed(K);
    rank_r1 = xmalloc((size_t)n_cap * sizeof(int));
    rank_r2 = xmalloc((size_t)n_cap * sizeof(int));
    rank_dst = xmalloc((size_t)n_cap * sizeof(int));
    memset(rank_r1, 0, (size_t)n_cap * sizeof(int));
    memset(rank_r2, 0, (size_t)n_cap * sizeof(int));
    memset(rank_dst, 0, (size_t)n_cap * sizeof(int));
    printf("===== 开销测量 =====\n");
    psi_print(stdout);
    printf("K=%d T=%d q=256 M=%d 丢包=%d%% 节点=源、中继1、中继2、目的\n", K, T, M, LOSS_PERCENT);
    printf("先发 %d 个 batch，再从少到多取前缀，找刚好能还原的个数。\n", n_cap);
    fflush(stdout);
    simulate_ex(src, K, T, code_seed, kLossSeed, kRecodeSeed, LOSS_PERCENT, n_cap, 0, -1,
                NULL, rank_r1, rank_r2, rank_dst, &dest);
    {
        int hist[M + 1];
        int r;
        memset(hist, 0, sizeof(hist));
        for (i = 0; i < n_cap; i++) {
            int rk = rank_dst[i];
            if (rk < 0) {
                rk = 0;
            }
            if (rk > M) {
                rk = M;
            }
            hist[rk]++;
        }
        printf("全部 %d 个发出 batch 的目的端秩个数（秩 0 到 %d）:\n", n_cap, M);
        for (r = 0; r <= M; r++) {
            printf("%d%c", hist[r], r == M ? '\n' : ' ');
            if (hist_out) {
                hist_out[r] = hist[r];
            }
        }
    }
    if (!try_prefix(src, K, T, code_seed, dest, n_cap, &bp, &inact)) {
        fprintf(stderr, "上界 %d 个 batch 仍不能还原\n", n_cap);
        goto done;
    }
    lo = 0;
    hi = n_cap;
    printf("搜索:");
    while (hi - lo > 1) {
        int mid = lo + (hi - lo) / 2;
        int step_bp = 0;
        int step_in = 0;
        int hit = try_prefix(src, K, T, code_seed, dest, mid, &step_bp, &step_in);
        printf(" %d%s", mid, hit ? "成" : "败");
        fflush(stdout);
        if (hit) {
            hi = mid;
            bp = step_bp;
            inact = step_in;
        } else {
            lo = mid;
        }
    }
    printf("\n");
    n_min = hi;
    if (!try_prefix(src, K, T, code_seed, dest, n_min, &bp, &inact)) {
        fprintf(stderr, "搜索得到的 %d 个 batch 不能还原\n", n_min);
        goto done;
    }
    if (n_min > 1 && try_prefix(src, K, T, code_seed, dest, n_min - 1, NULL, NULL)) {
        fprintf(stderr, "%d 个 batch 也能还原，最少个数不对\n", n_min - 1);
        goto done;
    }
    for (i = 0; i < n_min; i++) {
        rank_sum += rank_dst[i];
    }
    printf("刚好能还原的前 %d 个 batch 的秩分布（收到的系数矩阵，含整批丢失）\n", n_min);
    printf("  %-6s", "秩");
    for (i = 0; i <= M; i++) {
        printf(" %4d", i);
    }
    printf("   平均\n");
    print_rank_row("中继1", rank_r1, n_min);
    print_rank_row("中继2", rank_r2, n_min);
    print_rank_row("目的", rank_dst, n_min);
    printf("目的端秩总和 = %lld，平均秩 h_bar = %lld.%lld\n", rank_sum, rank_sum / n_min,
           (rank_sum * 10 / n_min) % 10);
    if (rank_sum > 0) {
        long long n_star_num = (long long)K * n_min;
        long long n_star = n_star_num / rank_sum;
        long long n_star_frac = (n_star_num * 10 / rank_sum) % 10;
        long long overhead = (rank_sum - (long long)K) * 100 / (long long)K;
        printf("理想批次数 K/h_bar = %lld.%lld\n", n_star, n_star_frac);
        printf("相对 K/h_bar 的开销 = %lld%%\n", overhead);
    }
    printf("BP 解出 %d 个源包，inactivation %d 个\n", bp, inact);
    ok = 1;
done:
    free_buckets(dest, n_cap);
    free(rank_r1);
    free(rank_r2);
    free(rank_dst);
    free(src);
    return ok;
}

int main(int argc, char **argv) {
    const char *in_path = "input.txt";
    const char *out_path = "recovered.txt";
    const char *psi_path = NULL;
    int optimize = 0;
    int iarg = 1;
    uint8_t toy[10 * 4];
    int i;
    Rng rng;
#ifdef _WIN32
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
#endif
    while (iarg < argc) {
        if (strcmp(argv[iarg], "--optimize") == 0) {
            optimize = 1;
            iarg++;
            continue;
        }
        if (strcmp(argv[iarg], "--psi") == 0) {
            if (iarg + 1 >= argc) {
                fprintf(stderr, "用法: bats_sim [--optimize | --psi 权重文件] [输入文件] [输出文件]\n");
                return 1;
            }
            psi_path = argv[iarg + 1];
            iarg += 2;
            continue;
        }
        if (argv[iarg][0] == '-') {
            fprintf(stderr, "用法: bats_sim [--optimize | --psi 权重文件] [输入文件] [输出文件]\n");
            return 1;
        }
        break;
    }
    if (optimize && psi_path) {
        fprintf(stderr, "--optimize 和 --psi 只能选一个\n");
        return 1;
    }
    if (argc - iarg > 2) {
        fprintf(stderr, "用法: bats_sim [--optimize | --psi 权重文件] [输入文件] [输出文件]\n");
        return 1;
    }
    gf_init();
    if (!gf_self_test() || !inactivation_self_test()) {
        return 1;
    }
    rng.s = 0x123456789ULL;
    for (i = 0; i < 10 * 4; i++) {
        toy[i] = rng_byte(&rng);
    }
    if (!run_silent(toy, 10, 4, 0xA11CEULL, 0, batches_needed(10)) ||
        !run_silent(toy, 10, 4, 0xA11CEULL, LOSS_PERCENT, batches_needed(10))) {
        fprintf(stderr, "线网自检失败\n");
        return 1;
    }
    printf("自检通过：GF(256)、BP 解不动后的 inactivation、无丢包线网、5%% 丢包线网\n");
    if (!run_example()) {
        return 1;
    }
    if (psi_path && !psi_load(psi_path)) {
        return 1;
    }
    {
        int hist[M + 1];
        int baseline_ok;
        int r;
        long long rank_n = 0;
        memset(hist, 0, sizeof(hist));
        baseline_ok = run_measure(hist);
        for (r = 0; r <= M; r++) {
            rank_n += hist[r];
        }
        if (optimize) {
            if (rank_n <= 0) {
                fprintf(stderr, "没有测到秩个数，无法优化\n");
                return 1;
            }
            if (!baseline_ok) {
                printf("当前分布没有在上界内还原，仍用上面的秩个数优化 Ψ\n");
            }
            printf("===== 按上面这次测量的秩分布优化 Ψ =====\n");
            if (!apply_psi(hist, K_MEASURE, "psi.txt")) {
                return 1;
            }
            printf("===== 用优化后的 Ψ 再测一次 =====\n");
            if (!run_measure(NULL)) {
                fprintf(stderr, "优化后的分布没有在上界内还原，改回写死的 Ψ\n");
                psi_reset();
                if (!run_measure(NULL)) {
                    return 1;
                }
            }
        } else if (!baseline_ok) {
            return 1;
        }
        if (iarg < argc) {
            in_path = argv[iarg++];
        }
        if (iarg < argc) {
            out_path = argv[iarg];
        }
        if (!run_file(in_path, out_path, optimize ? hist : NULL, optimize)) {
            return 1;
        }
    }
    return 0;
}
