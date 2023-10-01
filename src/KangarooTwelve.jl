module KangarooTwelve
# using SIMD

# See <https://www.fress.io/story/k12>
# For SIMD oppotunities, see <https://docs.rs/keccak/latest/src/keccak/lib.rs.html>.
# With native code, <https://lists.gnupg.org/pipermail/gcrypt-devel/2022-July/005360.html>
# and <https://github.com/XKCP/XKCP/blob/master/lib/low/KeccakP-1600/AVX2/KeccakP-1600-AVX2.s>
# for some hints as to what good assembly should look like.

import Base.Cartesian.@ntuple

## Keccak-p1600

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

function keccak_p1600(state::NTuple{25, UInt64}, ::Val{nrounds}=Val{12}()) where {nrounds}
    rol64(a, n) = (a << n) | (a >> (64 - n))
    # ~90ns (NB: the in-round times don't quite add up to this,
    # so they're not quite right, but they should still be a useful
    # indication)
    @inbounds for round in (25 - nrounds):24
        # θ (diffusion)
        C = @ntuple 5 i -> reduce(⊻, @ntuple 5 k -> state[i+5*(k-1)]) # ~12ns
        D = @ntuple 5 i -> C[mod1(i+4, 5)] ⊻ rol64(C[mod1(i+1, 5)], 1) # ~20ns
        state = @ntuple 25 i -> state[i] ⊻ D[mod1(i, 5)] # ~7ns
        # ρ (rotation) and π (lane permutation)
        state = @ntuple 25 i -> rol64(state[πs[i]], ρs[i]) # ~6ns
        # χ (intra-row bitwise combination, nonlinear)
        state = @ntuple 25 i -> state[i] ⊻ (~state[χs[i][1]] & state[χs[i][2]]) # ~10ns
        # ι (symmetry disruptor)
        state = Base.setindex(state, state[1] ⊻ ROUND_CONSTS_24[round], 1) # ~3ns
    end
    state
end

## TurboSHAKE

function ingest(state::NTuple{25, UInt64}, block::NTuple{rate, UInt64}) where {rate}
    state = @ntuple 25 i -> if i <= rate
        state[i] ⊻ block[i]
    else state[i] end
    keccak_p1600(state), rate
end

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

# 5% overhead compared to a UInt64 message
function ingest(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{U}) where {capacity, U<:Union{UInt32,UInt16,UInt8}}
    rate = 200 ÷ sizeof(U) - capacity ÷ (8 * sizeof(U))
    ratio = sizeof(UInt64)÷sizeof(U)
    for block in Iterators.partition(message, rate)
        state = if length(block) == rate
            @ntuple 25 i -> if i <= rate÷ratio
                state[i] ⊻ reduce(+, ntuple(
                    k -> UInt64(block[ratio*(i-1)+k]) << (8*sizeof(U)*(k-1)),
                    Val{ratio}()))
            else state[i] end
        else
            @ntuple 25 i ->
                state[i] ⊻ reduce(+, ntuple(
                    k -> UInt64(get(block, ratio*(i-1)+k, zero(U))) << (8*sizeof(U)*(k-1)),
                    Val{ratio}()))
        end |> keccak_p1600
    end
    state, fld1(mod1(length(message), rate), sizeof(U))
end

function pad(state::NTuple{25, UInt64}, ::Val{capacity}, finalblk::Int, delimsufix::UInt8) where {capacity}
    rate = 25 - capacity ÷ 64
    state = Base.setindex(state, state[finalblk] ⊻ UInt64(delimsufix), finalblk)
    state = Base.setindex(state, state[rate-1] ⊻ UInt64(0x80), rate-1)
    keccak_p1600(state)
end

function squeeze(::Type{NTuple{count, U}}, state::NTuple{25, UInt64}, ::Val{capacity}=Val{CAPACITY}()) where {capacity, count, U<:Unsigned}
    rate = 25 - capacity ÷ 64
    if count * sizeof(U) > rate * sizeof(UInt64)
        chunksize = rate * sizeof(UInt64) ÷ sizeof(U)
        chunkfirst = squeeze(NTuple{chunksize, U}, state, Val{capacity}())
        chunkrest = squeeze(NTuple{count - chunksize, U}, keccak_p1600(state), Val{capacity}())
        ntuple(i -> if i <= chunksize
                   chunkfirst[i]
               else
                   chunkrest[i - chunksize]
               end, Val{count}())
    elseif U == UInt128
        ntuple(i -> UInt128(state[2*(i-1)+1]) << 64 + state[2*i], Val{count}())
    elseif U == UInt64
        ntuple(i -> state[i], Val{count}())
    else
        byteratio = sizeof(UInt64)÷sizeof(U)
        ntuple(i -> (state[fld1(i, byteratio)] >> (8 * sizeof(U) * mod(-i, byteratio))) % U, Val{count}())
    end
end

squeeze(::Type{U}, state::NTuple{25, UInt64}, capacity::Val=Val{CAPACITY}()) where {U <: Unsigned} =
    squeeze(NTuple{1, U}, state, capacity) |> first

function squeeze!(output::Vector{UInt64}, state::NTuple{25, UInt64}, ::Val{capacity}=Val{CAPACITY}()) where {capacity}
    rate = 25 - capacity ÷ 64
    if length(output) <= rate
        output .= state[1:length(output)]
    else
        index = 1
        while index < length(output)
            index == 1 || (state = keccak_p1600(state))
            bsize = min(rate, length(output) - index + 1)
            output[index:index+bsize-1] .= state[1:bsize]
            index += bsize
        end
    end
    output
end

function turboshake(output::Type, # <:Unsigned or NTuple{n, <:Unsigned}
                    message::AbstractVector{<:Union{UInt64, UInt32, UInt16, UInt8}},
                    delimsufix::UInt8=0x80, capacity::Val = Val{CAPACITY}())
    state, finalblk = ingest(EMPTY_STATE, capacity, message)
    state = pad(state, capacity, finalblk, delimsufix)
    squeeze(output, state, capacity)
end

# ## KangarooTwelve

const BLOCK_SIZE = 8192
const CAPACITY = 256

const K12_SUFFIXES = (one=0x07, many=0x06, leaf=0x0b)

struct Trunk{rate}
    state::NTuple{25, UInt64}
    growth::NTuple{rate, UInt64}
    byte::Int # index
end

Trunk{rate}() where {rate} = Trunk(EMPTY_STATE, ntuple(_ -> zero(UInt64), Val{rate}()), 1)
Trunk() = Trunk{25 - CAPACITY ÷ 64}()

function overwrite(larger::Ubig, smaller::Usmall, byte::Int=1) where {Ubig <: Unsigned, Usmall <: Unsigned}
    # REVIEW should `htol` need to be used here?
    Ubig == Usmall && return smaller
    shift = 8 * (sizeof(Ubig) - byte - sizeof(Usmall) + 1)
    mask = ~(Ubig(typemax(Usmall)) << shift)
    value = Ubig(smaller) << shift
    larger & mask | value
end

function overwrite(larger::NTuple{size, Ubig}, smaller::Usmall, byte::Int=1) where {size, Ubig <: Unsigned, Usmall <: Unsigned}
    @boundscheck byte + sizeof(Usmall)-1 <= sizeof(Ubig) * size ||
        throw(BoundsError(larger, fld1(byte + sizeof(Usmall)-1, sizeof(Ubig))))
    largeindex = fld1(byte, sizeof(Ubig))
    if byte % sizeof(Ubig) == 1 || fld1(byte + sizeof(Usmall), sizeof(Ubig)) == largeindex
        # If `smaller` will fit neatly into the next `larger`, we can
        # just `overwrite` the targeted `larger` element with `smaller`.
        ntuple(i -> if i == largeindex
                          overwrite(larger[i], smaller, mod1(byte, sizeof(Ubig)))
                      else larger[i] end,
                      Val{size}())
    else
        # It doesn't fit, which means we'll need to split it up.
        # We can't actually use `overwrite(::Ubig, ::Usmall, byte)` as
        # we run into the potential of needing a UInt48 etc.
        shift = 8 * ((byte + sizeof(Usmall)-1) % sizeof(Ubig))
        mask1 = ~(Ubig(typemax(Usmall)) >> shift)
        value1 = Ubig(smaller) >> shift
        mask2 = ~(Ubig(typemax(Usmall)) << (8*sizeof(Ubig) - shift))
        value2 = Ubig(smaller) << (8*sizeof(Ubig) - shift)
        ntuple(i -> if i == largeindex
                   larger[i] & mask1 | value1
               elseif i == largeindex + 1
                   larger[i] & mask2 | value2
               else larger[i] end,
               Val{size}())
    end
end

function ingest(trunk::Trunk{rate}, entry::U) where {rate, U <: Union{UInt64, UInt32, UInt16, UInt8}}
    if trunk.byte <= 8 * rate - sizeof(U)+1 # Fits within `growth`
        growth = overwrite(trunk.growth, entry, trunk.byte)
        nextbyte = trunk.byte + sizeof(U)
        if nextbyte == 8 * rate + 1 # Rollover
            # TODO change `ntuple(_ -> zero(UInt64), Val{rate}())` to just
            # `growth`, once we've verified that this is behaving correctly.
            Trunk(ingest(trunk.state, growth) |> first,
                  ntuple(_ -> zero(UInt64), Val{rate}()), 1)
        else
            Trunk(trunk.state, growth, nextbyte)
        end
    else # Wrap-around needed, this also means we must be on the last byte of `growth`.
        growth = trunk.growth
        glast, gfirst = overwrite((growth[end], zero(UInt64)), entry,
                                  mod1(trunk.byte, sizeof(eltype(growth))))
        # See <https://github.com/JuliaLang/julia/issues/15276> for why the `let` is needed.
        growth = let growth=growth; ntuple(i -> if i == rate glast else growth[i] end, Val{rate}()) end
        state = ingest(trunk.state, growth) |> first
        # REVIEW if this seems like it might affect performance, we could
        # try changing `zero(UInt64)` to `growth[i]` so that less bytes need
        # to be modified. With `zero(UInt64)` it's easier to debug/inspect though.
        growth = ntuple(i -> if i == 1 gfirst else zero(UInt64) end, Val{rate}())
        Trunk(state, growth, mod1(trunk.byte + sizeof(U), sizeof(eltype(growth))))
    end
end

function k12_singlethreaded_allocfree(message::Vector{UInt64})
    if length(message) <= BLOCK_SIZE÷8
        return turboshake(UInt128, message, K12_SUFFIXES.one)
    end
    trunk = Trunk()
    rate = length(trunk.growth)
    # REVIEW is this the right way to initially ingest?
    for i in 1:rate
        trunk = ingest(trunk, message[i])
    end
    trunk = ingest(trunk, 0xc000000000000000)
    n, rest = 0, view(message, BLOCK_SIZE÷8+1:length(message))
    for block in Iterators.partition(rest, BLOCK_SIZE÷8)
        ingest(trunk, turboshake(UInt32, block, K12_SUFFIXES.leaf))
    end
    trunk = ingest(ingest(ingest(trunk, UInt32(n)), 0xffff), 0x01)
    # Do extra roll of trunk keccak if needed (likely)
    state = pad(trunk.state, Val{CAPACITY}(), mod1(n, rate), K12_SUFFIXES.many)
    squeeze(UInt128, state, Val{CAPACITY}())
end

function k12_singlethreaded(message::Vector{UInt64})
    if length(message) <= BLOCK_SIZE÷8
        turboshake(UInt128, message, K12_SUFFIXES.one)
    else
        nodestar = message[1:BLOCK_SIZE÷8]
        rest = view(message, BLOCK_SIZE÷8+1:length(message))
        push!(nodestar, 0xc000000000000000)
        n = 0
        for block in Iterators.partition(rest, BLOCK_SIZE÷8)
            push!(nodestar, turboshake(UInt64, block, K12_SUFFIXES.leaf))
            n += 1
        end
        push!(nodestar, UInt64(n), 0x000000000000ffff)
        turboshake(UInt128, nodestar, K12_SUFFIXES.many) |> first
    end
end

function k12_singlethreaded(io::IO)
    block = Vector{UInt8}(undef, BLOCK_SIZE)
    if readbytes!(io, block) < BLOCK_SIZE
        # small
        return zero(UInt128)
    end
    nodestar = Vector{UInt64}(undef, BLOCK_SIZE÷8)
    copyto!(nodestar, 1, reinterpret(UInt64, block), 1, BLOCK_SIZE÷8)
    push!(nodestar, 0xc000000000000000)
    size, n = BLOCK_SIZE, 0
    while (size = readbytes!(io, block)) == BLOCK_SIZE
        push!(nodestar, turboshake(UInt64, reinterpret(UInt64, block), K12_SUFFIXES.leaf))
        n += 1
    end
    if size < BLOCK_SIZE
        nfill = 8 - size % 8
        copyto!(block, size + 1, zeros(UInt8, nfill), 1, nfill)
        push!(nodestar, turboshake(UInt64,
                                   reinterpret(UInt64, view(block, 1:size+nfill)),
                                   K12_SUFFIXES.leaf))
        n += 1
    end
    push!(nodestar, UInt64(n), 0x000000000000ffff)
    turboshake(UInt128, nodestar, K12_SUFFIXES.many)
end

# TODO: Multithreaded implementation + heuristic
const k12 = k12_singlethreaded

include("throughput.jl")

end
