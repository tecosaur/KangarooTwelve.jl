const EMPTY_STATE = @ntuple 25 _ -> zero(UInt64)

const ROUND_CONSTS_24 =
    (0x0000000000000001, 0x0000000000008082, 0x800000000000808a, 0x8000000080008000,
     0x000000000000808b, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
     0x000000000000008a, 0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
     0x000000008000808b, 0x800000000000008b, 0x8000000000008089, 0x8000000000008003,
     0x8000000000008002, 0x8000000000000080, 0x000000000000800a, 0x800000008000000a,
     0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008)

const ρs = UInt64.(
    (0, 44, 43, 21, 14, 28, 20, 3, 45, 61, 1, 6, 25,
     8, 18, 27, 36, 10, 15, 56, 62, 55, 39, 41, 2))

const πs =
    (1, 7, 13, 19, 25, 4, 10, 11, 17, 23, 2, 8, 14,
     20, 21, 5, 6, 12, 18, 24, 3, 9, 15, 16, 22)

const χs =
    (( 2,  3), ( 3,  4), ( 4,  5), ( 5,  1), ( 1,  2), ( 7,  8),
     ( 8,  9), ( 9, 10), (10,  6), ( 6,  7), (12, 13), (13, 14),
     (14, 15), (15, 11), (11, 12), (17, 18), (18, 19), (19, 20),
     (20, 16), (16, 17), (22, 23), (23, 24), (24, 25), (25, 21),
     (21, 22))

"""
    keccak_p1600(state::NTuple{25, UInt64}, ::Val{nrounds}=Val{12}())

Apply the Keccak-*p*[`nrounds`, 1600] permutation to `state`. This is formally
defined in [the Keccak reference](https://keccak.team/files/Keccak-reference-3.0.pdf)
 and formalised in [FIPS 202](https://csrc.nist.gov/pubs/fips/202/final).

# Extended help

This is a variable-round permutation, with up to 24 rounds, where the last round
performed matches the final Keccak-*f* permutation round.

Each round of permutation consists of five steps, termed θ, ρ, π, χ, and ι.
These all operate on a 200-byte state, viewed in a number of different
configurations.

```text
             ┌─┬─┬─┬─┬─┐
            ┌─┬─┬─┬─┬─┐┤
           ┌─┬─┬─┬─┬─┐┤┤
          ┌─┬─┬─┬─┬─┐┤┤┤
         ┌─┬─┬─┬─┬─┐┤┤┤┤
        ┌─┬─┬─┬─┬─┐┤┤┤┤┘
       ┌─┬─┬─┬─┬─┐┤┤┤┤┘
      ┌─┬─┬─┬─┬─┐┤┤┤┤┘
      ├─┼─┼─┼─┼─┤┤┤┤┘
      ├─┼─┼─┼─┼─┤┤┤┘
      ├─┼─┼─┼─┼─┤┤┘        ┌─┐ bit
      ├─┼─┼─┼─┼─┤┘         └─┘         ┌─┐
      └─┴─┴─┴─┴─┘                     ┌─┐┘
         state                       ┌─┐┘
                                    ┌─┐┘
                   ┌─┐             ┌─┐┘
                   ├─┤ column     ┌─┐┘
     row           ├─┤           ┌─┐┘
 ┌─┬─┬─┬─┬─┐       ├─┤          ┌─┐┘  lane
 └─┴─┴─┴─┴─┘       ├─┤          └─┘
                   └─┘
```

- **θ step** Compute the parity of five columns, and xor-diffuse their parity into nearby columns.
- **ρ step** Bitwise-rotate each of the 25 lanes by a different triangular number.
- **π step** Permute each lane in a fixed pattern.
- **χ step** Intra-row bitwise combination. This provides the non-linearity.
- **ι step** The first lane is (xor-)mixed with a LFSR sequence across rounds.
  This serves to disrupt the symmetry of the scheme.
"""
@inline function keccak_p1600(state::NTuple{25, <:Union{UInt64, <:Vec{<:Any, UInt64}}}, ::Val{nrounds}=Val{12}()) where {nrounds}
    # Inlining `roll64` makes this faster with SIMD, but also non-deterministic 😢
    rol64(a, n) = (a << n) | (a >> (64 - n))
    @inbounds for round in (25 - nrounds):24
        # θ (diffusion)
        C = @ntuple 5 i -> reduce(⊻, @ntuple 5 k -> state[i+5*(k-1)])
        D = @ntuple 5 i -> C[mod1(i+4, 5)] ⊻ rol64(C[mod1(i+1, 5)], 1)
        state = @ntuple 25 i -> state[i] ⊻ D[mod1(i, 5)]
        # ρ (rotation) and π (lane permutation)
        state = @ntuple 25 i -> rol64(state[πs[i]], ρs[i]) # the `rol64` SIMD issue occurs here
        # χ (intra-row bitwise combination, nonlinear)
        state = @ntuple 25 i -> state[i] ⊻ (~state[χs[i][1]] & state[χs[i][2]])
        # ι (symmetry disruptor)
        state = setindex(state, state[1] ⊻ ROUND_CONSTS_24[round], 1)
    end
    state
end
