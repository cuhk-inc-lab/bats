# BATS

这个目录是一份 BATS 编解码库，外加一个只在本进程内存里跑的四节点仿真。同事要接到自己的程序里，用库即可，不必带上仿真。

库做三件事：`bats_encode`、`bats_recode`、`bats_decode`。没有预编码，没有套接字，也不规定包头。`code_seed`、`batch_id`、`K`、`T` 由调用方自己带到对端。

## 原理

源数据切成 `K` 个包，每包 `T` 字节。

每次编码处理一个 batch：

1. 从度数分布 Ψ 抽出度数 `d`，再均匀选出 `d` 个源包。
2. 这 `d` 个包乘上一个 `d×M` 的随机矩阵 `G`，得到 `M = 8` 个编码载荷。
3. `G`、被选中的源包下标都不放进包里。编码端和译码端用同一个 `code_seed` 加上这个 batch 的 `batch_id`，各自把 `G` 算出来。
4. 第 `j` 个发出包的系数是 `M×M` 单位阵的第 `j` 列。载荷已经是编码结果。

中继只混合**同一个** `batch_id` 里已经收到的包。没收齐也可以混合。不同 batch 不能放进同一次再编码。

目的端不再编码。先做置信传播：当前秩等于度数就解出来。解不动时再做 inactivation。收齐并解出全部 `K` 个源包才算成功。

两端必须使用同一份 Ψ。不调用 `psi_set` 时用库内写好的分布。

## 接到别的程序

```bash
make lib
```

得到 `libbats.a`。调用方只包含 `bats_codec.h`：

```bash
gcc -std=c11 -I/path/to/bats caller.c /path/to/bats/libbats.a -o caller
```

库由 `field.c`、`source.c`、`dest.c`、`codec.c` 组成。`link.c`、`relay.c`、`sim.c`、`main.c` 是仿真，不要链进去。程序启动后调用一次 `bats_init`。

一个 batch 写出 `BATS_M`（8）个包：

| 缓冲区 | 大小 | 排布 |
| --- | --- | --- |
| `coeff` | `8 * 8` | 第 `j` 个包的系数从 `coeff[j * 8]` 开始 |
| `payload` | `8 * T` | 第 `j` 个包的载荷从 `payload[j * T]` 开始 |

源数据是连续的 `K * T` 字节。`bats_encode` 成功返回 `0`，参数不合法返回 `-1`。

`bats_recode` 吃进同一 batch 的 `n_in` 个包，写出 8 个新包。`n_in == 0` 时不写包。同一跳上后面的 batch 继续传入同一个 `recode_state`。输入和输出缓冲区可以重叠。

`bats_decode` 接收扁平的 `BatsSymbol` 数组，按 `batch_id` 分组后再解。解出全部 `K` 个包返回 `1`，方程还不够返回 `0`，参数不合法返回 `-1`。每个 `payload` 指针指向 `T` 字节，只在这次调用返回前有效。

```c
#include "bats_codec.h"

enum { K = 128, T = 32 };

void example(const uint8_t src[K * T]) {
    uint8_t coeff[BATS_M * BATS_M];
    uint8_t payload[BATS_M * T];
    uint8_t dst[K * T];
    BatsSymbol symbols[BATS_M];
    uint64_t code_seed = 1;
    uint32_t batch_id = 0;
    int j;

    bats_init();
    bats_encode(src, K, T, code_seed, batch_id, coeff, payload);

    for (j = 0; j < BATS_M; j++) {
        symbols[j].batch_id = batch_id;
        memcpy(symbols[j].coeff, coeff + j * BATS_M, BATS_M);
        symbols[j].payload = payload + j * T;
    }
    /* 返回 1 表示 dst 里是原来的 K 个包。实际链路上通常要多个 batch。 */
    bats_decode(symbols, BATS_M, K, T, code_seed, dst, NULL, NULL);
}
```

下一个 batch 把 `batch_id` 加 1 再编一次。对端用同一个 `code_seed`。丢了哪些包由调用方决定，库不模拟丢包。

## 仿真

仿真是源 → 中继 1 → 中继 2 → 目的，包不出进程。丢包发生在入队之前。

```bash
make
./bats_sim input.txt recovered.txt
./bats_sim --optimize input.txt recovered.txt
./bats_sim --psi psi.txt input.txt recovered.txt
```

`--optimize` 和 `--psi` 只能选一个。还原结果与输入一致时返回 `0`，否则返回 `1`。程序还会对 `K = 128` 搜索刚好能还原的 batch 数，并打印每跳秩分布。

`psi_opt` 是同一套度数分布优化器，可以单独跑。

Windows 上如果有 `tools/tcc/tcc.exe`，用 `build.bat`。

## 文件

| 文件 | 作用 |
| --- | --- |
| `bats_codec.h` | 给外部程序的接口 |
| `codec.c` | `bats_encode`、`bats_recode`、`bats_decode` |
| `field.c` | GF(256) 和 batch 随机数 |
| `source.c` | 度数分布和 `G` |
| `dest.c` | 组方程、置信传播、inactivation |
| `link.c` `relay.c` `sim.c` `main.c` | 只给仿真用 |
| `psi_opt.c` | 按测量到的秩分布求 Ψ |
