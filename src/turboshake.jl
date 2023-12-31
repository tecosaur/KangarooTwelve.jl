# Utility functions

cap2rate(capacity::Integer) = 25 - capacity ÷ 64
rate2cap(rate::Integer) = (25 - rate) * 64

# Because `reinterpret(NTuple{...}, x)` and friends only works on 1.10+
@static if VERSION >= v"1.10-beta"
    ntupleinterpret(::Type{T}, x::U) where {T<:Unsigned, U<:Unsigned} =
        reinterpret(NTuple{sizeof(U) ÷ sizeof(T), T}, x)
    ntupleinterpret(::Type{T}, x::NTuple{N, U}) where {T<:Unsigned, N, U<:Unsigned} =
        reinterpret(NTuple{N * sizeof(U) ÷ sizeof(T), T}, x)
    ntupleinterpret(::Type{T}, x::NTuple{N, NTuple{M, T}}) where {T<:Unsigned, N, M} =
        reinterpret(NTuple{N * M, T}, x)
    uinterpret(::Type{T}, x::Union{U, NTuple{N, U}}) where {T<:Unsigned, N, U<:Unsigned} =
        reinterpret(T, x)
else
    @generated function ntupleinterpret(::Type{T}, x::U) where {T<:Unsigned, U<:Unsigned}
        Expr(:tuple, ntuple(i -> :(x >> $(8 * (i - 1) * sizeof(T)) % $T), sizeof(U) ÷ sizeof(T))...)
    end
    @generated function ntupleinterpret(::Type{T}, x::NTuple{N, U}) where {T<:Unsigned, N, U<:Unsigned}
        if sizeof(T) <= sizeof(U)
            Expr(:tuple, ntuple(i -> :(x[$(fld1(i, sizeof(U) ÷ sizeof(T)))] >> $(8 * ((i - 1) % (sizeof(U) ÷ sizeof(T))) * sizeof(T)) % $T),
                                N * sizeof(U) ÷ sizeof(T))...)
        else
            Expr(:tuple, ntuple(i ->
                Expr(:call, :+, ntuple(j -> :((x[$(sizeof(T) ÷ sizeof(U) * (i-1) + j)] % $T) << $(8 * sizeof(U) * (j-1))), sizeof(T) ÷ sizeof(U))...),
                                (N * sizeof(U)) ÷ sizeof(T))...)
        end
    end
    @generated function ntupleinterpret(::Type{T}, x::NTuple{N, NTuple{M, T}}) where {T<:Unsigned, N, M}
        Expr(:tuple, ntuple(i -> :(x[$(fld1(i, M))][$(mod1(i, M))]), N * M)...)
    end
    uinterpret(::Type{T}, x::Union{U, NTuple{N, U}}) where {T<:Unsigned, N, U<:Unsigned} =
        first(ntupleinterpret(T, x))
end

# Ingestion, padding, and squeezing

"""
    absorb(state::NTuple{25, UInt64}, block::NTuple{rate, UInt64})

Ingest a single `block` of input (with implied `rate`) into `state`.

The first `rate` elements of the state are `xor`'d `block`, and then the state
is permuted with `keccak_p1600`.
"""
function absorb(state::NTuple{25, UInt64}, block::NTuple{rate, UInt64}) where {rate}
    state = @ntuple 25 i -> if i <= rate
        state[i] ⊻ block[i]
    else state[i] end
    keccak_p1600(state)
end

"""
    absorb(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{<:Unsigned})

Ingest `message` into `state`, with the rate calculated based on `capacity`.

This breaks `message` into rate-sized blocks and then absorbs them (as per
`absorb(state, ::NTuple{rate, UInt64})`) in turn.
"""
function absorb(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{UInt64}) where {capacity}
    rate = cap2rate(capacity)
    loopend = length(message) - length(message) % rate
    for offset in 1:rate:loopend
        block = @view message[offset:offset-1+rate]
        state = keccak_p1600(@ntuple 25 i -> state[i] ⊻ if i <= rate block[i] else zero(UInt64) end)
    end
    block = @view message[loopend+1:end]
    return @ntuple 25 i -> state[i] ⊻ get(block, i, zero(UInt64))
end

# Unfortunately we do have to write a second SIMD-capable version since the
# block construction is a little different to the linear version.
function absorb(state::NTuple{25, Vec{N, UInt64}}, ::Val{capacity}, messages::NTuple{N, <:AbstractVector{UInt64}}) where {N, capacity}
    rate = cap2rate(capacity)
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

# ~5% overhead compared to a UInt64 message
function absorb(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{U}) where {capacity, U<:Union{UInt32,UInt16,UInt8}}
    rate = 200 ÷ sizeof(U) - capacity ÷ (8 * sizeof(U))
    ratio = sizeof(UInt64)÷sizeof(U)
    for block in partition(message, rate)
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

"""
    pad(state::NTuple{25, UInt64}, ::Val{capacity}, lastbyte::UInt, delimsuffix::UInt8)

Perform "padding" of the `lastbyte` byte of `state` with `delimsuffix`.

This is the final modification to `state` in [`turboshake`](@ref).
"""
function pad(state::NTuple{25, UInt64}, ::Val{capacity}, lastbyte::UInt, delimsuffix::UInt8) where {capacity}
    rate = cap2rate(capacity)
    last_uint64 = fld1(lastbyte, sizeof(UInt64))
    padbyte = state[last_uint64] ⊻ hton((delimsuffix % UInt64) << (8 * (7 - (lastbyte-1) % 8)))
    state = setindex(state, padbyte, last_uint64)
    state = setindex(state, state[rate] ⊻ (0x80 % UInt64) << 56, rate)
    keccak_p1600(state)
end

# For SIMD results
pad(states::NTuple{25, Vec{N, UInt64}}, capacity::Val, lastbyte::NTuple{N, UInt}, delimsuffix::UInt8) where {N} =
    ntuple(i -> pad(ntuple(k -> states[k][i], Val{25}()), capacity, lastbyte[i], delimsuffix), Val{N}())

"""
    squeeze(outtype::Type, state::NTuple{25, UInt64}, ::Val{capacity})

Squeeze an `outtype` out of `state`. `outtype` can be an `Unsigned` type or an
unsigned `NTuple`.
"""
function squeeze(::Type{NTuple{count, U}}, state::NTuple{25, UInt64}, ::Val{capacity}) where {capacity, count, U<:Unsigned}
    rate = cap2rate(capacity)
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
        ntuple(i -> state[2*(i-1)+1] + UInt128(state[2*i]) << 64, Val{count}())
    elseif U == UInt64
        ntuple(i -> state[i], Val{count}())
    else
        byteratio = sizeof(UInt64)÷sizeof(U)
        ntuple(i -> (state[fld1(i, byteratio)] >> (8 * sizeof(U) * mod(i-1, byteratio))) % U, Val{count}())
    end
end

squeeze(::Type{U}, state::NTuple{25, UInt64}, capacity::Val) where {U <: Unsigned} =
    squeeze(NTuple{1, U}, state, capacity) |> first

# For SIMD results
squeeze(U::Type, states::NTuple{N, NTuple{25, UInt64}}, capacity::Val) where {N} =
    ntuple(i -> squeeze(U, states[i], capacity), Val{N}())

# `squeeze!` isn't actually called anywhere, but it seems nice to have for completeness.
"""
    squeeze!(output::Vector{UInt64}, state::NTuple{25, UInt64}, ::Val{capacity})

Squeeze `state` into `output`.
"""
function squeeze!(output::AbstractVector{U}, state::NTuple{25, UInt64}, ::Val{capacity}) where {capacity, U<:Unsigned}
    rate = cap2rate(capacity)
    if length(output) * sizeof(U) <= rate * 8
        ssize = div(length(output) * sizeof(U), 8, RoundUp)
        oslice = ntupleinterpret(U, state[1:ssize])
        output .= oslice[1:length(output)]
    else
        index = 1
        while index < length(output)
            index == 1 || (state = keccak_p1600(state))
            osize = min(rate * 8 ÷ sizeof(U), length(output) - index + 1)
            ssize = div(osize * sizeof(U), 8, RoundUp)
            oslice = ntupleinterpret(U, state[1:ssize])
            output[index:index+osize-1] .= oslice[1:osize]
            index += osize
        end
    end
    output
end

# Turboshake frontend

# These functions are manually inlined in order to facilitate the compiler
# determining that no allocations are actually needed.

"""
    turboshake(output::Type, message::AbstractVector{<:UInt8to64},
               delimsuffix::UInt8=0x80, ::Val{capacity} = Val{$CAPACITY}())

Produce an `output` (`Unsigned` or `NTuple{n, <:Unsigned}`) value, by
performing TurboSHAKE on `message` with a certain `capacity` and `delimsuffix`.

See also: [`absorb`](@ref), [`pad`](@ref), and [`squeeze`](@ref).
"""
@inline function turboshake(output::Type, # <:Unsigned or NTuple{n, <:Unsigned}
                    message::AbstractVector{<:UInt8to64},
                    delimsuffix::UInt8=0x80, ::Val{capacity} = Val{CAPACITY}()) where {capacity}
    state = absorb(EMPTY_STATE, Val{capacity}(), message)
    # It might seem like `mod1` would make sense here, but for *whatever reason*
    # that seems to cause allocations. This took several hours to pinpoint.
    lastbyte = (length(message) * sizeof(eltype(message))) % UInt % (200 - capacity ÷ 8) + 1
    state = pad(state, Val{capacity}(), lastbyte, delimsuffix)
    squeeze(output, state, Val{capacity}())
end

# SIMD version
@inline function turboshake(output::Type, # <:Unsigned or NTuple{n, <:Unsigned}
                    message::NTuple{N, <:AbstractVector{<:UInt8to64}},
                    delimsuffix::UInt8=0x80, ::Val{capacity} = Val{CAPACITY}()) where {N, capacity}
    empty_state = ntuple(_ -> Vec(ntuple(_ -> zero(UInt64), Val{N}())), Val{25}())
    state = absorb(empty_state, Val{capacity}(), message)
    # This `lastindex` expression allocates. Why!?
    lastindex = ntuple(i -> (length(message[i]) * sizeof(eltype(message[i]))) % UInt % (200 - capacity ÷ 8) + 1, Val{N}())
    state = pad(state, Val{capacity}(), lastindex, delimsuffix)
    squeeze(output, state, Val{capacity}())
end
