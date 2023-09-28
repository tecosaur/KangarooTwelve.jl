module KangarooTwelve
# using SIMD

# See <https://www.fress.io/story/k12>
# For SIMD oppotunities, see <https://docs.rs/keccak/latest/src/keccak/lib.rs.html>.
# With native code, <https://lists.gnupg.org/pipermail/gcrypt-devel/2022-July/005360.html>
# and <https://github.com/XKCP/XKCP/blob/master/lib/low/KeccakP-1600/AVX2/KeccakP-1600-AVX2.s>
# for some hints as to what good assembly should look like.

import Base.Cartesian.@ntuple

## Keccak-p1600

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

const EMPTY_STATE = Tuple(zeros(UInt64, 25))

function keccak_p1600(state::NTuple{25, UInt64}, ::Val{nrounds}=Val{12}()) where {nrounds}
    rol64(a, n) = (a << n) | (a >> (64 - n))
    # ~90ns
    @inbounds for round in (25 - nrounds):24
        # θ (diffusion)
        C = @ntuple 5 i -> reduce(⊻, @ntuple 5 k -> state[i+5*(k-1)]) # ~12ns
        D = @ntuple 5 i -> C[mod1(i+4, 5)] ⊻ rol64(C[mod1(i+1, 5)], 1) # ~20ns
        state = @ntuple 25 i -> state[i] ⊻ D[mod1(i, 5)] # ~7ns
        # ρ (rotation) and π (lane permutation)
        state = @ntuple 25 i -> rol64(state[πs[i]], ρs[i]) # ~6ns
        # χ (intra-row bitwise combination, nonlinear)
        state = ( # ~10ns
            state[ 1] ⊻ (~state[ 2] & state[ 3]),
            state[ 2] ⊻ (~state[ 3] & state[ 4]),
            state[ 3] ⊻ (~state[ 4] & state[ 5]),
            state[ 4] ⊻ (~state[ 5] & state[ 1]),
            state[ 5] ⊻ (~state[ 1] & state[ 2]),
            state[ 6] ⊻ (~state[ 7] & state[ 8]),
            state[ 7] ⊻ (~state[ 8] & state[ 9]),
            state[ 8] ⊻ (~state[ 9] & state[10]),
            state[ 9] ⊻ (~state[10] & state[ 6]),
            state[10] ⊻ (~state[ 6] & state[ 7]),
            state[11] ⊻ (~state[12] & state[13]),
            state[12] ⊻ (~state[13] & state[14]),
            state[13] ⊻ (~state[14] & state[15]),
            state[14] ⊻ (~state[15] & state[11]),
            state[15] ⊻ (~state[11] & state[12]),
            state[16] ⊻ (~state[17] & state[18]),
            state[17] ⊻ (~state[18] & state[19]),
            state[18] ⊻ (~state[19] & state[20]),
            state[19] ⊻ (~state[20] & state[16]),
            state[20] ⊻ (~state[16] & state[17]),
            state[21] ⊻ (~state[22] & state[23]),
            state[22] ⊻ (~state[23] & state[24]),
            state[23] ⊻ (~state[24] & state[25]),
            state[24] ⊻ (~state[25] & state[21]),
            state[25] ⊻ (~state[21] & state[22]),
        )
        # ι (symmetry disruptor)
        state = Base.setindex(state, state[1] ⊻ ROUND_CONSTS_24[round], 1) # ~3ns
    end
    state
end

## TurboSHAKE

function ingest(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{UInt64}) where {capacity}
    rate = 25 - capacity ÷ 64
    for block in Iterators.partition(message, rate)
        state = if length(block) == rate
            @ntuple 25 i -> if i <= rate
                state[i] ⊻ block[i]
            else state[i] end
        else
            @ntuple 25 i -> state[i] ⊻ get(block, i, zero(UInt64))
        end |> keccak_p1600
    end
    state, mod1(length(message), rate)
end

for U in (UInt32, UInt16, UInt8)
    @eval function ingest(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{$U}) where {capacity}
        rate = 200 ÷ sizeof(U) - capacity ÷ (8 * sizeof(U))
        for block in Iterators.partition(message, rate)
            state = if length(block) == rate
                @ntuple 25 i -> if i <= rate÷$(sizeof(UInt64)÷sizeof(U))
                    state[i] ⊻ reduce(+, @ntuple $(sizeof(UInt64)÷sizeof(U)) k ->
                        UInt64(block[$(sizeof(UInt64)÷sizeof(U))*(i-1)+k]) <<
                            (8*sizeof(U)*(k-1)))
                else state[i] end
            else
                @ntuple 25 i ->
                    state[i] ⊻ reduce(+, @ntuple $(sizeof(UInt64)÷sizeof(U)) k ->
                        UInt64(get(block, $(sizeof(UInt64)÷sizeof(U))*(i-1)+k,
                                   zero(U))) << (8*sizeof(U)*(k-1)))
            end |> keccak_p1600
        end
        state, fld1(mod1(length(message), rate), sizeof(U))
    end
end

function pad(state::NTuple{25, UInt64}, ::Val{capacity}, finalblk::Int, delimsufix::UInt8) where {capacity}
    rate = 25 - capacity ÷ 64
    state = Base.setindex(state, state[finalblk] ⊻ UInt64(delimsufix), finalblk)
    state = Base.setindex(state, state[rate-1] ⊻ UInt64(0x80), rate-1)
    keccak_p1600(state)
end

function squeeze(state::NTuple{25, UInt64}, ::Val{capacity}, ::Val{output}) where {capacity, output}
    rate = 25 - capacity ÷ 64
    if output ÷ 64 <= rate
        state[1:output÷64]
    else
        squeeze!(zeros(UInt64, output÷64), state, Val{capacity}())
    end
end

function squeeze!(output::Vector{UInt64}, state::NTuple{25, UInt64}, ::Val{capacity}) where {capacity}
    rate = 25 - capacity ÷ 64
    if length(output) <= rate
        output .= state[1:length(output)]
    else
        index = 1
        while index < length(output)
            index == 1 || (state = keccak_p1600(state))
            bsize = min(rate, length(output) - index + 1)
            output[index:index+bsize] .= state[1:bsize]
            index += bsize
        end
    end
    output
end

function turboshake(message::AbstractVector{<:Union{UInt64, UInt32, UInt16, UInt8}},
                    delimsufix::UInt8=0x80, capacity::Val = Val{CAPACITY}(),
                    output::Val = (function (::Val{c}) where {c} Val{c÷2}() end)(capacity))
    state, finalblk = ingest(EMPTY_STATE, capacity, message)
    state = pad(state, capacity, finalblk, delimsufix)
    squeeze(state, capacity, output)
end

# ## KangarooTwelve

const BLOCK_SIZE = 8192
const CAPACITY = 256

function k12_singlethreaded(message::Vector{UInt64})
    high, low =  if length(message) <= BLOCK_SIZE÷8
        turboshake(message, 0x07)
    else
        nodestar = message[1:BLOCK_SIZE÷8]
        rest = view(message, BLOCK_SIZE÷8+1:length(message))
        push!(nodestar, 0xc000000000000000)
        n = 0
        for block in Iterators.partition(rest, BLOCK_SIZE÷8)
            push!(nodestar, turboshake(block, 0x0b, Val{CAPACITY}(),
                                       Val(64)) |> first)
            n += 1
        end
        push!(nodestar, UInt64(n), 0x000000000000ffff)
        turboshake(nodestar, 0x06, Val{CAPACITY}(), Val{128}())
    end
    UInt128(high) << 64 + low
end

function k12_singlethreaded(io::IO)
    block = Vector{UInt8}(undef, BLOCK_SIZE)
    if readbytes!(io, block) < BLOCK_SIZE
        # small
        return
    end
    nodestar = Vector{UInt64}(undef, BLOCK_SIZE÷8)
    copyto!(nodestar, 1, reinterpret(UInt64, block), 1, BLOCK_SIZE÷8)
    push!(nodestar, 0xc000000000000000)
    size, n = BLOCK_SIZE, 0
    while (size = readbytes!(io, block)) == BLOCK_SIZE
        push!(nodestar, turboshake(reinterpret(UInt64, block), 0x0b,
                                   Val{CAPACITY}(), Val{64}()) |> first)
        n += 1
    end
    if size < BLOCK_SIZE
        nfill = 8 - size % 8
        copyto!(block, size + 1, zeros(UInt8, nfill), 1, nfill)
        push!(nodestar, turboshake(
            reinterpret(UInt64, view(block, 1:size+nfill)),
            0x0b, Val{CAPACITY}(), Val{64}()) |> first)
        n += 1
    end
    push!(nodestar, UInt64(n), 0x000000000000ffff)
    high, low = turboshake(nodestar, 0x06)
    UInt128(high) << 64 + low
end

# TODO: Multithreaded implementation + heuristic
const k12 = k12_singlethreaded

## Testing

function throughput(::typeof(keccak_p1600), size) # 2 = 1KiB, 2 = 32KiB, 4 = 1MiB...
    state = Tuple(zeros(UInt64, 25))
    rounds = round(Int, 32^size / 200)
    start = time()
    for _ in 1:rounds
        state = keccak_p1600(state)
    end
    stop = time()
    push!([], state) # prevent the call from being optimised away
    println(" Keccak-p[1600,12] throughput: ~$(round(Int, 32^size / (stop - start) / 1024^2)) MiB/s")
end

function throughput(::typeof(turboshake), size)
    message = rand(UInt64, round(Int, 32^size) ÷ 8)
    start = time()
    x = turboshake(message)
    stop = time()
    push!([], x) # prevent the call from being optimised away
    println(" TurboSHAKE128 throughput: ~$(round(Int, 32^size / (stop - start) / 1024^2)) MiB/s")
end

function throughput(::typeof(k12_singlethreaded), size)
    message = rand(UInt64, round(Int, 32^size) ÷ 8)
    start = time()
    x = k12_singlethreaded(message)
    stop = time()
    push!([], x) # prevent the call from being optimised away
    println(" K12 (singlethreaded) throughput: ~$(round(Int, 32^size / (stop - start) / 1024^2)) MiB/s")
end

function throughput(step::Symbol, size)
    func = if step == :keccak
        keccak_p1600
    elseif step == :turboshake
        turboshake
    elseif step == :k12_singlethread
        k12_singlethreaded
    elseif step == :k12
        k12
    else
        error("I don't know that step of the k12 algorithm.")
    end
    throughput(func, size)
end

function throughputs(size=6)
    for step in [:keccak, :turboshake, :k12_singlethread]
        throughput(step, size)
    end
end

end
