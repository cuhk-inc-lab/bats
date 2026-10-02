# BATS line-network simulator

In-process simulator for one flow: source → relay 1 → relay 2 → destination.
Packets never leave the process. Loss is applied before a packet is queued.
There is no precoding, no socket, and no port.

Each batch draws a degree `d` from a fixed distribution, selects `d` source
packets, and emits `M = 8` coded packets. `G` is not carried in the packet.
Relays mix only packets that share a `batch_id`. The destination runs belief
propagation, then inactivation.

## Build

```bash
make
./bats_sim input.txt recovered.txt
```

On Windows, if `tools/tcc/tcc.exe` is present:

```bat
build.bat
bats_sim.exe input.txt recovered.txt
```

The program returns 0 when `recovered.txt` matches the input, and 1 otherwise.
It also searches for the smallest batch count that recovers `K = 128` and
prints the per-hop rank histogram.

By default the degree distribution stays the built-in one. `--optimize` keeps
the destination rank counts from the `K = 128` measurement even when that
run cannot recover, solves (P1) with maximum degree 128, and measures again.
File recovery then solves another distribution capped at that file's own
packet count, and leaves `psi.txt` as the `K = 128` result. `--psi` uses a
weight file instead. Pass one of those options, not both.

```bash
./bats_sim input.txt recovered.txt
./bats_sim --optimize input.txt recovered.txt
./bats_sim --psi psi.txt input.txt recovered.txt
```

`psi_opt` is the same optimizer as a separate program, if you already have
rank counts and want to write a weight file yourself.

## Calling from another program

Include `bats_codec.h` only. Link `field.c`, `source.c`, `dest.c`, and `codec.c`
(`make lib` builds `libbats.a`). Those four files are the codec. The line
network stays in `link.c`, `relay.c`, `sim.c`, and `main.c`. Call `bats_init` once.

`bats_encode` writes one batch: `BATS_M` coefficient vectors and `BATS_M`
payloads, each payload `T` bytes. `bats_recode` mixes packets that already
share a `batch_id`; pass the same `recode_state` for later batches on that
hop. `bats_decode` takes a flat list of `BatsSymbol` and groups them by
`batch_id`. A return value of 1 means all `K` source packets were solved.

```c
uint8_t coeff[BATS_M * BATS_M];
uint8_t payload[BATS_M * 32];
uint64_t recode_state = 1;
int n_out = 0;

bats_init();
bats_encode(src, K, 32, code_seed, batch_id, coeff, payload);
bats_recode(coeff, payload, BATS_M, 32, &recode_state, coeff, payload, &n_out);
```

## Files

| File | Role |
| --- | --- |
| `source.c` | Degree distribution and `G` |
| `relay.c` | Relay: recode one batch. Relay 1 and relay 2 share this code |
| `dest.c` | Build batch equations, belief propagation, then inactivation |
| `sim.c` | Send each batch, then source → relay 1 → relay 2 → destination |
| `link.c` | In-memory links. Loss happens before enqueue |
| `field.c` | GF(256) and the batch RNG |
| `codec.c` | `bats_encode`, `bats_recode`, `bats_decode` |
| `bats_codec.h` | Header for callers outside this simulator |
| `main.c` | Self-check, small example, overhead search, file recovery |
