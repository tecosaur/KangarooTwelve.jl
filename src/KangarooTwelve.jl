module KangarooTwelve
# using SIMD

# See <https://www.fress.io/story/k12>
# For SIMD oppotunities, see <https://docs.rs/keccak/latest/src/keccak/lib.rs.html>.
# With native code, <https://lists.gnupg.org/pipermail/gcrypt-devel/2022-July/005360.html>
# and <https://github.com/XKCP/XKCP/blob/master/lib/low/KeccakP-1600/AVX2/KeccakP-1600-AVX2.s>
# for some hints as to what good assembly should look like.

import Base.Cartesian.@ntuple

const ROUND_CONSTS_24 =
    [0x0000000000000001, 0x0000000000008082, 0x800000000000808a, 0x8000000080008000,
     0x000000000000808b, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
     0x000000000000008a, 0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
     0x000000008000808b, 0x800000000000008b, 0x8000000000008089, 0x8000000000008003,
     0x8000000000008002, 0x8000000000000080, 0x000000000000800a, 0x800000008000000a,
     0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008]

const ρs = UInt64[
    1,  3,  6,  10, 15, 21, 28, 36, 45, 55,  2, 14,
    27, 41, 56,  8, 25, 43, 62, 18, 39, 61, 20, 44]

const πs = Int[
    11,  8, 12, 18, 19, 4,  6,  17, 9, 22, 25, 5,
    16, 24, 20, 14, 13, 3, 21, 15, 23, 10,  7, 2]

# The `keccak_p1600!` function currently takes ~290ns / 1500 cycles
# to execute. The K12 paper states about their C implementation:
#   "On an Intel Core i5-6500 (Skylake), we measured that one evaluation of
#    Keccak-p[1600, nr = 12] takes about 450 cycles"
# This suggests that we /should/ be able to get this down to ~90ns per
# evaluation.

function keccak_p1600!(state::Vector{UInt64}, ::Val{nrounds}=Val{12}()) where {nrounds}
    @boundscheck length(state) == 25 || throw(ArgumentError("`state` must be of length exactly 25 (not $(length(state)))"))
    rol64(a, n) = (a >> (64 - (n % 64))) | (a << (n % 64))
    @inbounds for rc in view(ROUND_CONSTS_24, (25-nrounds):24)
        # θ (diffusion)
        C = @ntuple 5 i -> reduce(⊻, @ntuple 5 k -> state[i+5*(k-1)]) # ~1ns
        D = @ntuple 5 i -> C[mod1(i+4, 5)] ⊻ rol64(C[mod1(i+1, 5)], 1) # ~0ns
        for i in 1:25 # ~120ns
            state[i] ⊻= D[mod1.(i, 5)]
        end
        # # ρ (rotation) and π (lane permutation)
        last = state[πs[end]]
        for i in 1:24 # ~100ns
            j = state[πs[i]]
            state[πs[i]] = rol64(last, ρs[i])
            last = j
        end
        # # χ (intra-row bitwise combination, nonlinear)
        for j in 0:5:20 # ~50ns
            C = @ntuple 5 i -> state[j + i]
            for i in 1:5
                state[j + i] = C[i] ⊻ (~C[mod1(i+1, 5)] & C[mod1(i+2, 5)])
            end
        end
        # ι (symmetry disruptor)
        state[1] ⊻= rc # ~1ns
    end
    state
end

end
