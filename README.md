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

## Files

| File | Role |
| --- | --- |
| `source.c` | Source: sample a degree, build `G`, emit `M` packets |
| `relay.c` | Relay: recode one batch. Relay 1 and relay 2 share this code |
| `dest.c` | Destination: belief propagation, then inactivation |
| `sim.c` | Source → relay 1 → relay 2 → destination |
| `link.c` | In-memory links. Loss happens before enqueue |
| `field.c` | GF(256) and the batch RNG |
| `main.c` | Self-check, small example, overhead search, file recovery |
