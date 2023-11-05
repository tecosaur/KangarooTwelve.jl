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

A Keccak state that keeps track of the last *lane* updated.

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

"""
    absorb(sponge::AbstractSponge, block::AbstractVector{<:Unsigned})

Ingest each element of `block` into `sponge`.
"""
function absorb(sponge::AbstractSponge, block::AbstractVector{U}) where {U <: UInt8to64}
    # REVIEW optimize? This is just a quick hack
    for b in block
        sponge = absorb(sponge, b)
    end
    sponge
end

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

ByteSponge{rate}() where {rate} = ByteSponge{rate}(EMPTY_STATE, 1)

Base.convert(::Type{ByteSponge}, (; state, lane)::Sponge{rate}) where {rate} =
    ByteSponge{rate}(state, sizeof(UInt64) * (lane - 1) + 1)

function Base.convert(::Type{Sponge}, (; state, byte)::ByteSponge{rate}) where {rate}
    if (byte - 1) % sizeof(UInt64) == 0
        Sponge{rate}(state, (byte - 1) ÷ sizeof(UInt64) + 1)
    else
        (@noinline () -> throw(ArgumentError("ByteSponge could not be cleanly converted to a Sponge as the current byte did not cleanly align with a lane")))()
    end
end

pad((; state, byte)::ByteSponge{rate}, delimsuffix::UInt8) where {rate} =
    ByteSponge{rate}(pad(state, Val{rate2cap(rate)}(), byte, delimsuffix), 1)

"""
    subxor(larger::Ubig, smaller::Usmall, byte::UInt=1)

xor `larger` with `smaller`, lining the start of `smaller` up with `byte` of
`larger`.
"""
function subxor(larger::Ubig, smaller::Usmall, byte::UInt=1) where {Ubig <: Unsigned, Usmall <: Unsigned}
    Ubig == Usmall && return larger ⊻ smaller
    shift = 8 * (byte - 1)
    value = (smaller % Ubig) << shift
    larger ⊻ value
end

"""
    absorb(sponge::ByteSponge, x::Union{<:Unsigned, NTuple{N, <:Unsigned}})

Ingest `x` into `sponge`.
"""
function absorb((; state, byte)::ByteSponge{rate}, x::UInt8) where {rate}
    lane, lbyte = (byte-1) ÷ sizeof(UInt64) + 1, (byte-1) % sizeof(UInt64) + 1
    state = setindex(state, subxor(state[lane], x, lbyte), lane)
    byte += 1
    if lane == rate && lbyte == sizeof(UInt64)
        state, byte = keccak_p1600(state), 1
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
