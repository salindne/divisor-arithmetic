# smalljac wall-clock comparison

[`bench.c`](bench.c) times smalljac's genus-2 imaginary/ramified group law
(`hecurve_g2_compose` = add, `hecurve_g2_square` = double) on the affine path
(`ctx = NULL`, one field inversion per op), so it can be put next to this
crate's `g2::ramified::not_char2` `cargo bench` numbers. See the
"Wall-clock comparison with smalljac" section of the top-level
[`README.md`](../../README.md) for the results and the caveats.

smalljac and ff_poly are **not vendored** here (they are GPL; this crate is
MIT). Download them yourself:

- smalljac v4.1.3 and ff_poly v1.2.7: <https://math.mit.edu/~drew/smalljac.html>
- GMP (e.g. `brew install gmp` on macOS)

## Building on x86-64 Linux

ff_poly's hand-written x86-64 assembly works as-is there:

```sh
# build libff_poly.a per its Makefile, then:
cc -O3 -std=gnu99 -DSMALLJAC_GENUS=2 \
   -I/path/to/ff_poly -I/path/to/smalljac bench.c \
   hecurve.c hecurve1.c hecurve2.c mpzutil.c prime.c \
   /path/to/libff_poly.a -lgmp -o sjbench
./sjbench
```

## Building on arm64 (Apple Silicon)

ff_poly and smalljac use x86-64 inline asm (`mulq`, `bsrq`, `bsfq`, …). Replace
it with portable `__uint128_t` / `__builtin_*` equivalents (semantics identical;
clang lowers the 128-bit ops to native `mul`/`umulh`/`clz`/`ctz` on aarch64).

**1. Replace `ff_poly/asm.h`** with these portable macros:

```c
#define _asm_div_q_q(q,r,x,y)        do { __uint128_t _n=((__uint128_t)(unsigned long)(r)<<64)|(unsigned long)(x); (q)=(unsigned long)(_n/(unsigned long)(y)); (r)=(unsigned long)(_n%(unsigned long)(y)); } while(0)
#define _asm_mult_1_1(z1,z0,x0,y0)   do { __uint128_t _p=(__uint128_t)(unsigned long)(x0)*(unsigned long)(y0); (z0)=(unsigned long)_p; (z1)=(unsigned long)(_p>>64); } while(0)
#define _asm_mult_2_2_1(z1,z0,x1,x0,y0) do { __uint128_t _p=(__uint128_t)(unsigned long)(x0)*(unsigned long)(y0); (z0)=(unsigned long)_p; (z1)=(unsigned long)(_p>>64)+(unsigned long)(x1)*(unsigned long)(y0); } while(0)
#define _asm_addto_2_2(z1,z0,x1,x0)  do { __uint128_t _s=(((__uint128_t)(unsigned long)(z1)<<64)|(unsigned long)(z0))+(((__uint128_t)(unsigned long)(x1)<<64)|(unsigned long)(x0)); (z0)=(unsigned long)_s; (z1)=(unsigned long)(_s>>64); } while(0)
#define _asm_addto_2_1(z1,z0,x0)     do { __uint128_t _s=(((__uint128_t)(unsigned long)(z1)<<64)|(unsigned long)(z0))+(unsigned long)(x0); (z0)=(unsigned long)_s; (z1)=(unsigned long)(_s>>64); } while(0)
#define _asm_addto_3_3(z2,z1,z0,x2,x1,x0) do { __uint128_t _s=(__uint128_t)(unsigned long)(z0)+(unsigned long)(x0); (z0)=(unsigned long)_s; _s=(__uint128_t)(unsigned long)(z1)+(unsigned long)(x1)+(unsigned long)(_s>>64); (z1)=(unsigned long)_s; (z2)=(unsigned long)(z2)+(unsigned long)(x2)+(unsigned long)(_s>>64); } while(0)
#define _asm_addto_3_2(z2,z1,z0,x1,x0)    do { __uint128_t _s=(__uint128_t)(unsigned long)(z0)+(unsigned long)(x0); (z0)=(unsigned long)_s; _s=(__uint128_t)(unsigned long)(z1)+(unsigned long)(x1)+(unsigned long)(_s>>64); (z1)=(unsigned long)_s; (z2)=(unsigned long)(z2)+(unsigned long)(_s>>64); } while(0)
#define _asm_subfrom_2_2(z1,z0,x1,x0) do { __uint128_t _d=(((__uint128_t)(unsigned long)(z1)<<64)|(unsigned long)(z0))-(((__uint128_t)(unsigned long)(x1)<<64)|(unsigned long)(x0)); (z0)=(unsigned long)_d; (z1)=(unsigned long)(_d>>64); } while(0)
#define _asm_inc_2(z1,z0)            do { if(++(z0)==0UL) (z1)++; } while(0)
```

(plus the `_asm_mult_3_2_1` / `_asm_mult_3_2_2` / `_asm_square_3_2` helpers,
which are already written in terms of the macros above — keep them as-is.)

**2. In both `cstd.h` files** (the ff_poly copy *and* the smalljac copy) replace
the two bit-scan helpers:

```c
static inline unsigned long _asm_highbit (unsigned long x) { return 63UL - __builtin_clzl(x); }
static inline unsigned long _asm_lowbit  (unsigned long x) { return __builtin_ctzl(x); }
```

Then build with the same `cc` line as above. ff_poly is then limited to
`FF_BITS = 57`, which is why the matched-width comparison uses a 56-bit prime.
