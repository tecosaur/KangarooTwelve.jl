module KangarooTwelve
using Base.Threads
using SIMD
using Mmap

# See <https://www.fress.io/story/k12>
# For SIMD oppotunities, see <https://docs.rs/keccak/latest/src/keccak/lib.rs.html>.
# With native code, <https://lists.gnupg.org/pipermail/gcrypt-devel/2022-July/005360.html>
# and <https://github.com/XKCP/XKCP/blob/master/lib/low/KeccakP-1600/AVX2/KeccakP-1600-AVX2.s>
# for some hints as to what good assembly should look like.
#
# The performance should be similar to <https://keccak.team/sw_performance.html>

export k12

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

"""
    keccak_p1600(state::NTuple{25, UInt64}, ::Val{nrounds}=Val{12}())

Apply the Keccak-p[`nrounds`, 1600] permutation to `state`. This is formally
defined in [the Keccak reference](https://keccak.team/files/Keccak-reference-3.0.pdf)
 and formalised in [FIPS 202](https://csrc.nist.gov/pubs/fips/202/final).
"""
function keccak_p1600(state::NTuple{25, <:Union{UInt64, <:Vec{<:Any, UInt64}}}, ::Val{nrounds}=Val{12}()) where {nrounds}
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
        state = Base.setindex(state, state[1] ⊻ ROUND_CONSTS_24[round], 1)
    end
    state
end

## TurboSHAKE

const CAPACITY = 256

"""
    ingest(state::NTuple{25, UInt64}, block::NTuple{rate, UInt64})

Ingest a single `block` of input (with implied `rate`) into `state` or `trunk`.

The first `rate` elements of the state are `xor`'d `block`, and then the state
is permuted with `keccak_p1600`.
"""
function ingest(state::NTuple{25, UInt64}, block::NTuple{rate, UInt64}) where {rate}
    state = @ntuple 25 i -> if i <= rate
        state[i] ⊻ block[i]
    else state[i] end
    keccak_p1600(state)
end

"""
    ingest(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{<:Unsigned})

Ingest `message` into `state`, with the rate calculated based on `capacity`.

This breaks `message` into rate-sized blocks and then ingests them (as per
`ingest(state, ::NTuple{rate, UInt64})`) in turn.
"""
function ingest(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{UInt64}) where {capacity}
    rate = 25 - capacity ÷ 64
    for block in Iterators.partition(message, rate)
        state = if length(block) == rate
            @ntuple 25 i -> if i <= rate
                state[i] ⊻ block[i]
            else state[i] end
        else
            return @ntuple 25 i -> state[i] ⊻ get(block, i, zero(UInt64))
        end |> keccak_p1600
    end
    state
end

# Unfortunately we do have to write a second SIMD-capable version since the
# block construction is a little different to the linear version.
function ingest(state::NTuple{25, Vec{N, UInt64}}, ::Val{capacity}, messages::NTuple{N, <:AbstractVector{UInt64}}) where {N, capacity}
    rate = 25 - capacity ÷ 64
    msglengths = map(length, messages)
    for pos in 1:rate:maximum(msglengths)
        state = if all(>=(pos+rate), msglengths)
            @ntuple 25 i -> if i <= rate
                state[i] ⊻ Vec{N, UInt64}(ntuple(k -> messages[k][pos-1+i], Val{N}()))
            else state[i] end
        else
            return @ntuple 25 i -> if i <= rate
                state[i] ⊻ Vec{N, UInt64}(ntuple(k -> get(messages[k], pos-1+i, zero(UInt64)), Val{N}()))
            else state[i] end
        end |> keccak_p1600
    end
    state
end

# 5% overhead compared to a UInt64 message
function ingest(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{U}) where {capacity, U<:Union{UInt32,UInt16,UInt8}}
    rate = 200 ÷ sizeof(U) - capacity ÷ (8 * sizeof(U))
    ratio = sizeof(UInt64)÷sizeof(U)
    for block in Iterators.partition(message, rate)
        state = if length(block) == rate
            @ntuple 25 i -> if i <= rate÷ratio
                state[i] ⊻ reduce(+, ntuple(
                    k -> (block[ratio*(i-1)+k] % UInt64) << (8*sizeof(U)*(k-1)),
                    Val{ratio}()))
            else state[i] end
        else
            return @ntuple 25 i ->
                state[i] ⊻ reduce(+, ntuple(
                    k -> (get(block, ratio*(i-1)+k, zero(U)) % UInt64) << (8*sizeof(U)*(k-1)),
                    Val{ratio}()))
        end |> keccak_p1600
    end
    state
end

function pad(state::NTuple{25, UInt64}, ::Val{capacity}, lastbyte::Int, delimsuffix::UInt8) where {capacity}
    rate = 25 - capacity ÷ 64
    last_uint64 = fld1(lastbyte, sizeof(UInt64))
    padbyte = state[last_uint64] ⊻ hton((delimsuffix % UInt64) << (8 * (7 - (lastbyte-1) % 8)))
    state = Base.setindex(state, padbyte, last_uint64)
    state = Base.setindex(state, state[rate] ⊻ (0x80 % UInt64) << 56, rate)
    keccak_p1600(state)
end

# For SIMD results
pad(states::NTuple{25, Vec{N, UInt64}}, capacity, lastbyte::NTuple{N, Int64}, delimsuffix::UInt8) where {N} =
    ntuple(i -> pad(ntuple(k -> states[k][i], Val{25}()), capacity, lastbyte[i], delimsuffix), Val{N}())

"""
    squeeze(outtype::Type, state::NTuple{25, UInt64}, ::Val{capacity}=Val{CAPACITY}())

Squeeze an `outtype` out of `state`. `outtype` can be an `Unsigned` type or an
unsigned `NTuple`.
"""
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

# For SIMD results
squeeze(U::Type, states::NTuple{N, NTuple{25, UInt64}}, capacity::Val=Val{CAPACITY}()) where {N} =
    ntuple(i -> squeeze(U, states[i], capacity), Val{N}())
squeeze(U::Type, states::NTuple{25, Vec{N, UInt64}}, capacity::Val=Val{CAPACITY}()) where {N} =
    squeeze(U, unwrap_simd(states), capacity)

unwrap_simd(vals::NTuple{V, Vec{N, T}}) where {V, N, T} =
    ntuple(n -> ntuple(v -> vals[v][n], Val{V}()), Val{N}())

# `squeeze!` isn't actually called anywhere, but it seems nice to have for completeness.
"""
    squeeze!(output::Vector{UInt64}, state::NTuple{25, UInt64}, ::Val{capacity}=Val{CAPACITY}())

Squeeze `state` into `output`.
"""
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
                    delimsuffix::UInt8=0x80, ::Val{capacity} = Val{CAPACITY}()) where {capacity}
    state = ingest(EMPTY_STATE, Val{capacity}(), message)
    # It might seem like `mod1` would make sense here, but for *whatever reason*
    # that seems to cause allocations. This took several hours to pinpoint.
    lastbyte = (length(message) * sizeof(eltype(message))) % (200 - capacity ÷ 8) + 1
    state = pad(state, Val{capacity}(), lastbyte, delimsuffix)
    squeeze(output, state, Val{capacity}())
end

# SIMD version
function turboshake(output::Type, # <:Unsigned or NTuple{n, <:Unsigned}
                    message::NTuple{N, <:AbstractVector{<:Union{UInt64, UInt32, UInt16, UInt8}}},
                    delimsuffix::UInt8=0x80, ::Val{capacity} = Val{CAPACITY}()) where {N, capacity}
    empty_state = ntuple(_ -> Vec(ntuple(_ -> zero(UInt64), Val{N}())), Val{25}())
    state = ingest(empty_state, Val{capacity}(), message)
    # This `lastindex` expression allocates. Why!?
    lastindex = ntuple(i -> (length(message[i]) * sizeof(eltype(message[i]))) % (200 - capacity ÷ 8) + 1, Val{N}())
    state = pad(state, Val{capacity}(), lastindex, delimsuffix)
    squeeze(output, state, Val{capacity}())
end

## Trunk

struct Trunk{rate}
    state::NTuple{25, UInt64}
    byte::Int # index
end

Trunk{rate}() where {rate} = Trunk{rate}(EMPTY_STATE, 1)
Trunk() = Trunk{25 - CAPACITY ÷ 64}()

function Trunk{rate}(zeroblock::AbstractVector{U}) where {rate, U<:Unsigned}
    capacity = (25 - rate) * 64
    trunk = if length(zeroblock) * sizeof(U) == BLOCK_SIZE
        zerostate = ingest(EMPTY_STATE, Val{capacity}(), reinterpret(UInt64, zeroblock))
        Trunk{rate}(zerostate, 1 + BLOCK_SIZE % (8 * rate))
    else
        zerostate = ingest(EMPTY_STATE, Val{capacity}(), zeroblock)
        Trunk{rate}(zerostate, 1 + (length(zeroblock) * sizeof(U)) % (8 * rate))
    end
    ingest(trunk, K12_ZEROBLOCK_SUFFIX)
end

Trunk(zeroblock::AbstractVector{<:Unsigned}) =
    Trunk{25 - CAPACITY ÷ 64}(zeroblock)

ingest(trunk::Trunk{rate}, block::NTuple{rate, UInt64}) where {rate} =
    Trunk{rate}(ingest(trunk.state, block), 1)

function pad(trunk::Trunk{rate}, delimsuffix::UInt8) where {rate}
    (; state, byte) = trunk
    if byte != 1
        state = keccak_p1600(state)
    end
    Trunk{rate}(pad(state, Val{(25 - rate) * 64}(), byte, delimsuffix), 1)
end

squeeze(T::Type, trunk::Trunk{rate}) where {rate} =
    squeeze(T, trunk.state, Val{(25 - rate) * 64}())

## KangarooTwelve

const BLOCK_SIZE = 8192

const K12_SUFFIXES = (one=0x07, many=0x06, leaf=0x0b)

"""
    subxor(larger::Ubig, smaller::Usmall, byte::Int=1)
    subxor(larger::NTuple{size, Ubig}, smaller::Usmall, byte::Int=1)

xor `larger` with `smaller`, lining the start of `smaller` up with `byte` of
`larger`. When `larger` is a `NTuple`, `small` will be split across entries if
needed.
"""
function subxor(larger::Ubig, smaller::Usmall, byte::Int=1) where {Ubig <: Unsigned, Usmall <: Unsigned}
    Ubig == Usmall && return larger ⊻ smaller
    shift = 8 * (byte - 1)
    value = (smaller % Ubig) << shift
    larger ⊻ value
end

function subxor(larger::NTuple{size, Ubig}, smaller::Usmall, byte::Int=1) where {size, Ubig <: Unsigned, Usmall <: Unsigned}
    @boundscheck byte + sizeof(Usmall)-1 <= sizeof(Ubig) * size ||
        throw(BoundsError(larger, fld1(byte + sizeof(Usmall)-1, sizeof(Ubig))))
    largeindex = fld1(byte, sizeof(Ubig))
    if byte % sizeof(Ubig) == 1 || fld1(byte + sizeof(Usmall)-1, sizeof(Ubig)) == largeindex
        Base.setindex(larger,
                      subxor(larger[largeindex], smaller, mod1(byte, sizeof(Ubig))),
                      largeindex)
    else
        # It doesn't fit, which means we'll need to split it up.
        # We can't actually use `subxor(::Ubig, ::Usmall, byte)` as
        # we run into the potential of needing a UInt48 etc.
        shift = 8 * ((byte - 1) % sizeof(Ubig))
        overhang = 8 * (sizeof(Usmall) - (byte + sizeof(Usmall)-1) % sizeof(Ubig))
        value1 = (smaller % Ubig) << shift
        value2 = (smaller >> overhang) % Ubig
        larger = Base.setindex(larger, larger[largeindex] ⊻ value1, largeindex)
        Base.setindex(larger, larger[largeindex+1] ⊻ value2, largeindex+1)
    end
end

"""
    ingest(trunk::Trunk, entry::Unsigned)

Ingest a single `entry` into `trunk`. This updates the partial-block stored in
`trunk`, and when full merges it with the `state` and runs a `keccak_p1600`
round when full.
"""
function ingest(trunk::Trunk{rate}, entry::U) where {rate, U <: Union{UInt64, UInt32, UInt16, UInt8}}
    if trunk.byte <= 8 * rate - sizeof(U)+1 # Fits within `rate`
        state = subxor(trunk.state, entry, trunk.byte)
        nextbyte = trunk.byte + sizeof(U)
        if nextbyte == 8 * rate + 1 # Rollover
            Trunk{rate}(keccak_p1600(state), 1)
        else
            Trunk{rate}(state, nextbyte)
        end
    else # Wrap-around needed
        state = trunk.state
        slast, sfirst = subxor(
            (state[rate], zero(UInt64)), entry, mod1(trunk.byte, sizeof(UInt64)))
        state = Base.setindex(state, slast, rate)
        state = keccak_p1600(state)
        state = Base.setindex(state, sfirst, 1)
        Trunk{rate}(state, mod1(trunk.byte + sizeof(U), sizeof(UInt64)))
    end
end

"""
    ingest(trunk::Trunk, block::AbstractVector{<:Unsigned})

Ingest `block` into `trunk`. This is done by applying `turboshake` to `block`,
extracting a `UInt32`, and then ingesting that `UInt32` into `trunk`.
"""
function ingest(trunk::Trunk, block::AbstractVector{U}) where {U <: Union{UInt64, UInt32, UInt16, UInt8}}
    ingest(trunk, if U == UInt64
        turboshake(UInt32, block, K12_SUFFIXES.leaf)
    elseif length(block) == BLOCK_SIZE ÷ sizeof(U)
        turboshake(UInt32, reinterpret(UInt64, block), K12_SUFFIXES.leaf)
    else
        turboshake(UInt32, block, K12_SUFFIXES.leaf)
    end)
end

function ingest_length(trunk::Trunk, val::Int, ::Val{bufsize}=Val{8}()) where {bufsize}
    buffer = ntuple(_ -> 0x00, Val{bufsize}())
    point = 0
    while (val > 0)
        buffer = Base.setindex(buffer, UInt8(val % 2^8), point+=1)
        val ÷= 2^8
    end
    for i in point:-1:1
        trunk = ingest(trunk, buffer[i])
    end
    ingest(trunk, UInt8(point))
end

ingest(trunk::Trunk, x, xs...) = ingest(ingest(trunk, x), xs...)

const K12_ZEROBLOCK_SUFFIX = 0xc000000000000000
const SIMD_FACTOR = 4

function k12_singlethreaded_nosimd(message::AbstractVector{U}) where {U <: Union{UInt64, UInt32, UInt16, UInt8}}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return turboshake(UInt128, message, K12_SUFFIXES.one)
    end
    # Process the S₀ block
    rate = 25 - CAPACITY ÷ 64
    state0 = ingest(EMPTY_STATE, Val(CAPACITY), reinterpret(UInt64, view(message, 1:BLOCK_SIZE÷sizeof(U))))
    trunk = Trunk{rate}(state0, 1 + BLOCK_SIZE % (8 * rate))
    trunk = ingest(trunk, K12_ZEROBLOCK_SUFFIX)
    n, rest = 0, view(message, BLOCK_SIZE÷sizeof(U)+1:length(message))
    for block in Iterators.partition(rest, BLOCK_SIZE÷sizeof(U))
        out32 = turboshake(UInt32, block, K12_SUFFIXES.leaf)
        trunk, n = ingest(trunk, out32), n+1
    end
    trunk = ingest(ingest_length(trunk, n), 0xffff, 0x01)
    state = trunk.state
    # Do extra roll of trunk keccak if needed (likely)
    if trunk.byte != 1
        state = keccak_p1600(trunk.state)
    end
    state = pad(state, Val{CAPACITY}(), mod1(n, rate), K12_SUFFIXES.many)
    squeeze(UInt128, state, Val{CAPACITY}())
end

function k12_singlethreaded(message::AbstractVector{U}) where {U <: Union{UInt64, UInt32, UInt16, UInt8}}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return turboshake(UInt128, message, K12_SUFFIXES.one)
    end
    # Process the S₀ block
    rate = 25 - CAPACITY ÷ 64
    state0 = ingest(EMPTY_STATE, Val(CAPACITY), reinterpret(UInt64, view(message, 1:BLOCK_SIZE÷sizeof(U))))
    trunk = Trunk{rate}(state0, 1 + BLOCK_SIZE % (8 * rate))
    trunk = ingest(trunk, K12_ZEROBLOCK_SUFFIX)
    simd_block_size = SIMD_FACTOR * BLOCK_SIZE÷sizeof(U)
    message_simd_end = ((length(message) - BLOCK_SIZE÷sizeof(U)) ÷ simd_block_size) * simd_block_size + BLOCK_SIZE÷sizeof(U)
    n, centre = 0, view(message, BLOCK_SIZE÷sizeof(U)+1:message_simd_end)
    tail = view(message, message_simd_end+1:length(message))
    for block in Iterators.partition(centre, simd_block_size)
        block_u64 = reinterpret(UInt64, block)
        blocks = ntuple(i -> view(
            block_u64, ((i-1) * BLOCK_SIZE÷sizeof(UInt64) + 1):(i * BLOCK_SIZE÷sizeof(UInt64))),
                        SIMD_FACTOR)
        out_4u32 = turboshake(UInt32, blocks, K12_SUFFIXES.leaf)
        out1u64, out2u64 = reinterpret(NTuple{2, UInt64}, out_4u32)
        trunk, n = ingest(trunk, out1u64, out2u64), n + length(blocks)
    end
    for block in Iterators.partition(tail, BLOCK_SIZE÷sizeof(U))
        trunk, n = ingest(trunk, block), n+1
    end
    trunk = ingest(ingest_length(trunk, n), 0xffff, 0x01)
    state = trunk.state
    # Do extra roll of trunk keccak if needed (likely)
    if trunk.byte != 1
        state = keccak_p1600(trunk.state)
    end
    state = pad(state, Val{CAPACITY}(), mod1(n, rate), K12_SUFFIXES.many)
    squeeze(UInt128, state, Val{CAPACITY}())
end

# TODO: Better multithreaded implementation + heuristic
function k12_multithreaded(message::AbstractVector{U}) where {U <: Union{UInt64, UInt32, UInt16, UInt8}}
    chunks = Iterators.partition(message, 4 * BLOCK_SIZE÷sizeof(U))
    tasks = map(chunks) do chunk
        @spawn turboshake(UInt32, chunk, K12_SUFFIXES.leaf)
    end
    map(fetch, tasks)
end

function k12_singlethreaded(io::IO)
    block = Vector{UInt8}(undef, BLOCK_SIZE)
    if (size = readbytes!(io, block)) < BLOCK_SIZE
        return k12_singlethreaded(view(block, 1:size))
    end
    rate = 25 - CAPACITY ÷ 64
    trunk = Trunk{rate}()
    for i in 1:BLOCK_SIZE÷sizeof(UInt8)
        trunk = ingest(trunk, block[i])
    end
    trunk = ingest(trunk, K12_ZEROBLOCK_SUFFIX)
    n = 0
    bigblock = Vector{UInt8}(undef, 4 * BLOCK_SIZE)
    while readbytes!(io, block) > 0
        trunk, n = ingest(trunk, block), n+1
    end
    trunk = ingest(ingest_length(trunk, n), 0xffff, 0x01)
    state = trunk.state
    # Do extra roll of trunk keccak if needed (likely)
    if trunk.byte != 1
        state = keccak_p1600(trunk.state)
    end
    state = pad(state, Val{CAPACITY}(), mod1(n, rate), K12_SUFFIXES.many)
    squeeze(UInt128, state, Val{CAPACITY}())
end

"""
    k12(data::Union{IO, String, AbstractVector{Unsigned}}) -> UInt128

Hash `data` with the KangarooTwelve scheme.

This scheme presents a good balance of *Simplicity*, *Security*, and *Speed*.

# Extended help

The KangarooTwelve hashing scheme works by splitting the input data into ``n``
$BLOCK_SIZE-byte blocks (``S₀``, ``S₁``, …, ``Sₙ₋₁``) which are individually
hashed with [`TurboSHAKE128`](@ref turboshake) to produce "chaining
values" (CVs), which are put together and ingested to produce the final state.

```text
               ╭────╮ ╭────╮   ╭────╮ ╭────╮
               │ S₁ │ │ S₂ │   │Sₙ₋₂│ │Sₙ₋₁│
               ╰─┬──╯ ╰─┬──╯   ╰─┬──╯ ╰─┬──╯
                 │110   │110     │110   │110
                 ▼      ▼        ▼      ▼
╭────────╮110⁶²╭─┴──╮ ╭─┴──╮   ╭─┴──╮ ╭─┴──╮
│   S₀   ├─────┤ CV ├─┤ CV ├╴╍╶┤ CV ├─┤ CV ├─(n-1)(FFFF)(01)──▶─┤HASH│
╰────────╯     ╰────╯ ╰────╯   ╰────╯ ╰────╯
```

This scheme has been described as "leaves on a vine". The hashing of blocks
``S₁`` to ``Sₙ₋₁`` is embarassingly parallel, and can be accelerated with both
SIMD and multithreading. To ensure sufficient data to produce 128 bits of
entropy from a small number of blocks (each producing a `UInt32` CV), ``S₀`` is
included in full.
"""
function k12(data::AbstractVector{<:Unsigned})
    if true # length(data) < heuristic
        k12_singlethreaded(data)
    else
        k12_multithreaded(data)
    end
end

function k12(data::IO)
    try
        buf = Mmap.mmap(data)
        k12(data)
    catch _
        k12_multithreaded(data)
    end
end

k12(s::String) = k12(codeunits(s))
k12(s::AbstractString) = k12(String(s))

include("throughput.jl")

end
