# Most of the time, we'll be dealing with UInt64s,
# so it makes sense to specialise `Sponge` on that.

"""
    AbstractSponge

A sponge is a generalisation of hash functions and stream ciphers. It is in
essance an iterated construction of a variable-input variable-output function
based on a fixed length permutation and an internal state.

The way the sponge operates can be illustrated well with a diagram following the
state.

```text
 Message (padded)
  └───┬──────┬──────┬──────┐      ┊ ┌──────┬──────┬───▶ Output
      │      │      │      │      ┊ │      │      │
 ┌─┐  │  ╭─╮ │  ╭─╮ │  ╭─╮ │  ╭─╮ ┆ │  ╭─╮ │  ╭─╮ │
 │ │  ▼  │ │ ▼  │ │ ▼  │ │ ▼  │ │ ┆ │  │ │ │  │ │ │
r│0├──⨁─▶│ ├─⨁─▶│ ├─⨁─▶│ ├─⨁─▶│ ├─┼─┴─▶│ ├─┴─▶│ ├─┴─▶ ╍
 │ │     │f│    │f│    │f│    │f│ ┆    │f│    │f│
 ├─┤     │ │    │ │    │ │    │ │ ┆    │ │    │ │
c│0├────▶│ ├───▶│ ├───▶│ ├───▶│ ├─┼───▶│ ├───▶│ ├───▶ ╍
 └─┘     ╰─╯    ╰─╯    ╰─╯    ╰─╯ ┆    ╰─╯    ╰─╯
                        Absorbing ┆ Squeezing
```

First, the input is padded with a reversible padding rule and the state is
initialised to zero.
- In the absorbing phase, the input is processed in blocks of `r` bits, and xor'd with
  the first `r` bits of the state (leaving the remaining `c` bits unchanged), with
  absorptions interlieved with applications of the function `f`.
- In the squeezing phase, the first `r` bits of the state are returned as the output.
  Shoult more output be wanted, `f` can simply be applied to the state again.

While this construction is defined for any fixed-length permutation `f`, however
we specialise on [`keccak_p1600`](@ref).

It is possible to either incrementally either one lane at a time, or one byte at
a time. To handle these two cases we have the [`Sponge`](@ref) and [`ByteSponge`](@ref) subtypes.
"""
abstract type AbstractSponge{rate} end

"""
    Sponge{rate}

A Keccak state that keeps track of the current *lane*.

For more information on the sponge construction see [`AbstractSponge`](@ref).

See also: [`absorb`](@ref), [`pad`](@ref), and [`squeeze`](@ref).
"""
struct Sponge{rate} <: AbstractSponge{rate}
    state::NTuple{25, UInt64}
    lane::UInt
end

Sponge{rate}() where {rate} = Sponge{rate}(EMPTY_STATE, 1)
Sponge() = Sponge{RATE}()

function absorb((; state, lane)::Sponge{rate}, newlane::UInt64) where {rate}
    state = setindex(state, state[lane] ⊻ newlane, lane)
    if lane == rate
        Sponge{rate}(keccak_p1600(state), 1)
    else
        Sponge{rate}(state, lane + 1)
    end
end

function absorb((; state, lane)::Sponge{rate}, lanes::NTuple{N, UInt64}) where {rate, N}
    if rate == N && lane == 1
        state, lane = absorb(state, lanes), 1
    elseif lane + N <= rate
        state = @ntuple 25 i -> if lane <= i < lane + N
            state[i] ⊻ lanes[i - lane + 1]
        else state[i] end
        lane = lane + N
    else # Only works for N <= rate, but that should be fine
        state = @ntuple 25 i -> if lane <= i <= rate
            state[i] ⊻ lanes[i % UInt - lane + 1]
        else state[i] end
        state = keccak_p1600(state)
        state = @ntuple 25 i -> if i < N - rate + lane
            state[i] ⊻ lanes[(rate % UInt - lane) + i + 1]
        else state[i] end
        lane = (lane + N) % rate
    end
    Sponge{rate}(state, lane)
end

function absorb((; state, lane)::Sponge{rate}, lanes::AbstractVector{UInt64}) where {rate}
    if isempty(lanes)
        Sponge{rate}(state, lane)
    elseif lane == 1
        loopend = length(lanes) - length(lanes) % rate
        for offset in 1:rate:loopend
            block = @view lanes[offset:offset-1+rate]
            state = keccak_p1600(@ntuple 25 i ->
                state[i] ⊻ if i <= rate block[i] else zero(UInt64) end)
        end
        finalblock = @view lanes[loopend+1:end]
        state = @ntuple 25 i -> state[i] ⊻ get(finalblock, i, zero(UInt64))
        Sponge{rate}(state, (lane - 1 + length(lanes)) % rate + 1)
    else
        state = @ntuple 25 i -> if lane <= i <= min(rate, lane - 1 + length(lanes))
            state[i] ⊻ lanes[i - lane + 1]
        else state[i] end
        if length(lanes) >= rate - lane + 1
            state = keccak_p1600(state)
        end
        absorb(Sponge{rate}(state, min(rate, lane - 1 + length(lanes)) % rate + 1),
               view(lanes, (rate - lane + 2):lastindex(lanes)))
    end
end

"""
    squeeze(T::Type, sponge::AbstractSponge)

Squeeze a `T` out of `sponge`.
"""
squeeze(T::Type, (; state)::AbstractSponge{rate}) where {rate} =
    squeeze(T, state, Val{rate2cap(rate)}())

"""
    squeeze!(output::Vector{UInt64}, sponge::AbstractSponge)

Squeeze `sponge` into `output`.
"""
squeeze!(output::AbstractVector{<:Unsigned}, (; state)::AbstractSponge{rate}) where {rate} =
    squeeze!(output, state, Val{rate2cap(rate)}())

"""
    absorb(sponge::AbstractSponge, as::Type{<:Unsigned}, leaf::AbstractVector{<:Unsigned})

Ingest `leaf` into `sponge` by transforming it to an `as` via `turboshake`, and
absorbing that result.
"""
@inline absorb(sponge::AbstractSponge{rate}, T::Type, leaf::AbstractVector{U}) where {rate, U<:UInt8to64} =
    absorb(sponge, turboshake(T, leaf, K12_SUFFIXES.leaf))

# absorb(sponge::Sponge, x::Unsigned, xs::Unsigned...) = absorb(absorb(sponge, x), xs...)

# At the very end of KangarooTwelve, we want to absorb
# sub-UInt64 values. For this more limited capability,
# we need a second sponge type that can absorb single
# bytes of data.

"""
    ByteSponge{rate}

A Keccak state that keeps track of the last *byte* updated.

For more information on the sponge construction see [`AbstractSponge`](@ref).

See also: [`absorb`](@ref), [`pad`](@ref), and [`squeeze`](@ref).
"""
struct ByteSponge{rate} <: AbstractSponge{rate}
    state::NTuple{25, UInt64}
    byte::UInt
end

ByteSponge{rate}() where {rate} = ByteSponge{rate}(EMPTY_STATE, 0)

Base.convert(::Type{ByteSponge}, (; state, lane)::Sponge{rate}) where {rate} =
    ByteSponge{rate}(state, sizeof(UInt64) * (lane - 1))

function Base.convert(::Type{Sponge}, (; state, byte)::ByteSponge{rate}) where {rate}
    if byte % sizeof(UInt64) == 0
        Sponge{rate}(state, byte ÷ sizeof(UInt64) + 1)
    else
        (@noinline () -> throw(ArgumentError("ByteSponge could not be cleanly converted to a Sponge as the current byte did not cleanly align with a lane")))()
    end
end

pad((; state, byte)::ByteSponge{rate}, delimsuffix::UInt8) where {rate} =
    ByteSponge{rate}(pad(state, Val{rate2cap(rate)}(), byte+1, delimsuffix), 0)

"""
    subxor(larger::Ubig, smaller::Usmall, offset::UInt=0)

xor `larger` with `smaller`, lining the start of `smaller` up `offset` bytes
from the start of `larger`.
"""
function subxor(larger::Ubig, smaller::Usmall, offset::UInt=0) where {Ubig <: Unsigned, Usmall <: Unsigned}
    Ubig == Usmall && return larger ⊻ smaller
    shift = 8 * offset
    value = (smaller % Ubig) << shift
    larger ⊻ value
end

"""
    absorb(sponge::ByteSponge, x::Union{<:Unsigned, NTuple{N, <:Unsigned}})

Ingest `x` into `sponge`.
"""
function absorb((; state, byte)::ByteSponge{rate}, x::UInt8) where {rate}
    lane, lbyte = byte ÷ sizeof(UInt64) + 1, byte % sizeof(UInt64)
    state = setindex(state, subxor(state[lane], x, lbyte), lane)
    byte += 1
    if byte == rate * sizeof(UInt64)
        state, byte = keccak_p1600(state), 0
    end
    ByteSponge{rate}(state, byte)
end

function absorb(sponge::ByteSponge{rate}, x::U) where {rate, U<:Unsigned}
    for byte in ntupleinterpret(UInt8, x)
        sponge = absorb(sponge, byte)
    end
    sponge
end

function absorb(sponge::ByteSponge, x::NTuple{N, U}) where {N, U<:Unsigned}
    for val in x
        sponge = absorb(sponge, val)
    end
    sponge
end

"""
    absorb(sponge::AbstractSponge, block::AbstractVector{<:Unsigned})

Ingest each element of `block` into `sponge`.
"""
function absorb(sponge::ByteSponge{rate}, block::AbstractVector{U}) where {rate, U <: UInt8to64}
    if length(block) * sizeof(U) >= 8 * rate - sponge.byte >= 8
        start = min(firstindex(block), lastindex(block))
        aligngap = sponge.byte % sizeof(UInt64)
        if aligngap != 0
            for b in view(block, start:start+aligngap)
                sponge = absorb(sponge, b)
            end
            start += aligngap
        end
        remaining_bytes = (length(block) - aligngap) * sizeof(U)
        remaining_u64s = (remaining_bytes - remaining_bytes % sizeof(UInt64)) ÷ sizeof(U)
        sponge64 = convert(Sponge, sponge)
        sponge64 = absorb(sponge64, reinterpret(UInt64, view(block, start:start-1+remaining_u64s)))
        sponge = convert(ByteSponge, sponge64)
        for b in view(block, start+remaining_u64s:lastindex(block))
            sponge = absorb(sponge, b)
        end
    else
        for b in block
            sponge = absorb(sponge, b)
        end
    end
    sponge
end
