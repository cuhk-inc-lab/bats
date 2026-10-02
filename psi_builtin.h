#ifndef PSI_BUILTIN_H
#define PSI_BUILTIN_H

/* 写死的 Ψ。度数 1..8，权重和 256，平均度数 1186/256。
   source.c 和 psi_opt.c 都从这里取，避免两处各写一份。 */
#define PSI_BUILTIN_COUNT 8
#define PSI_BUILTIN_SUM 256
#define PSI_BUILTIN_W_INIT {12, 30, 40, 44, 40, 36, 30, 24}

#endif
