# Most of the time, we'll be dealing with UInt64s,
# so it makes sense to specialise `Sponge` on that.

abstract type AbstractSponge{rate} end

"""
    Sponge{rate}

A Keccak state that keeps track of the last *lane* updated.

See also: `ingest`, `pad`, and `squeeze`.
"""
struct Sponge{rate} <: AbstractSponge{rate}
    state::NTuple{25, UInt64}
    lane::UInt
end

Sponge{rate}() where {rate} = Sponge{rate}(EMPTY_STATE, 1)
Sponge() = Sponge{RATE}()

function ingest((; state, lane)::Sponge{rate}, newlane::UInt64) where {rate}
    state = setindex(state, state[lane] ⊻ newlane, lane)
    if lane == rate
        Sponge{rate}(keccak_p1600(state), 1)
    else
        Sponge{rate}(state, lane + 1)
    end
end

function ingest((; state, lane)::Sponge{rate}, lanes::NTuple{N, UInt64}) where {rate, N}
    if rate == N && lane == 1
        state, lane = ingest(state, lanes), 1
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

pad((; state, lane)::Sponge{rate}, delimsuffix::UInt8) where {rate} =
    Sponge{rate}(pad(state, Val{rate2cap(rate)}(), 8 * (lane-1) + 1, delimsuffix), 1)

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
    ingest(sponge::AbstractSponge, as::Type{<:Unsigned}, leaf::AbstractVector{<:Unsigned})

Ingest `leaf` into `sponge` by transforming it to an `as` via `turboshake`, and
ingesting that result.
"""
ingest(sponge::AbstractSponge{rate}, T::Type, leaf::AbstractVector{U}) where {rate, U<:UInt8to64} =
    ingest(sponge, turboshake(T, leaf, K12_SUFFIXES.leaf))

"""
    ingest(sponge::AbstractSponge, block::AbstractVector{<:Unsigned})

Ingest each element of `block` into `sponge`.
"""
function ingest(sponge::AbstractSponge, block::AbstractVector{U}) where {U <: UInt8to64}
    # REVIEW optimize? This is just a quick hack
    for b in block
        sponge = ingest(sponge, b)
    end
    sponge
end

# ingest(sponge::Sponge, x::Unsigned, xs::Unsigned...) = ingest(ingest(sponge, x), xs...)

# At the very end of KangarooTwelve, we want to ingest
# sub-UInt64 values. For this more limited capability,
# we need a second sponge type that can absorb single
# bytes of data.

"""
    ByteSponge{rate}

A Keccak state that keeps track of the last *byte* updated.

See also: `ingest`, `pad`, and `squeeze`.
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
    ingest(sponge::ByteSponge, x::Union{<:Unsigned, NTuple{N, <:Unsigned}})

Ingest `x` into `sponge`.
"""
function ingest((; state, byte)::ByteSponge{rate}, x::UInt8) where {rate}
    lane, lbyte = (byte-1) ÷ sizeof(UInt64) + 1, (byte-1) % sizeof(UInt64) + 1
    state = setindex(state, subxor(state[lane], x, lbyte), lane)
    byte += 1
    if lane == rate && lbyte == sizeof(UInt64)
        state, byte = keccak_p1600(state), 1
    end
    ByteSponge{rate}(state, byte)
end

function ingest(sponge::ByteSponge{rate}, x::U) where {rate, U<:Unsigned}
    for byte in ntupleinterpret(UInt8, x)
        sponge = ingest(sponge, byte)
    end
    sponge
end

function ingest(sponge::ByteSponge, x::NTuple{N, U}) where {N, U<:Unsigned}
    for val in x
        sponge = ingest(sponge, val)
    end
    sponge
end
