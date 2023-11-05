abstract type AbstractVine{rate} end

"""
    GerminatingVine{rate}

A `vine` that hasn't started to grow leaves yet (see the `k12` docstring).

See also: `Vine`, `ingest`, `ingest_length`, `finalise`.
"""
struct GerminatingVine{rate} <: AbstractVine{rate}
    trunk::ByteSponge{rate}
    leaf::ByteSponge{rate}
    nbytes::UInt
end

GerminatingVine{rate}() where {rate} =
    GerminatingVine{rate}(ByteSponge{rate}(), ByteSponge{rate}(), 0)

"""
    Vine{rate}

A `vine` with leaves (see the `k12` docstring).

See also: `GerminatingVine`, `ingest`, `ingest_length`, `finalise`.
"""
struct Vine{rate} <: AbstractVine{rate}
    trunk::Sponge{rate}
    leaf::ByteSponge{rate}
    nbytes::UInt
end

"""
    ingest(vine::AbstractVine, leaflet::AbstractVector{<:Unsigned})
    ingest(vine::AbstractVine, x::Unsigned)

Ingest `leaflet`/`x` into `vine`, this may go into the leaves or the trunk depending
on the vine type, and may cause the current leaf to be folded into the trunk,
and a new leaf grown.
"""
function ingest((; trunk, leaf, nbytes)::GerminatingVine{rate}, leaflet::AbstractVector{U}) where {rate, U<:Unsigned}
    leafletbytes = length(leaflet) * sizeof(U)
    if nbytes == 0 && leafletbytes == BLOCK_SIZE
        zerostate = ingest(EMPTY_STATE, Val{rate2cap(rate)}(), reinterpret(UInt64, leaflet))
        Vine{rate}(ingest(Sponge{rate}(zerostate, 1 + (BLOCK_SIZE ÷ 8) % rate), K12_ZEROBLOCK_SUFFIX),
                   leaf, BLOCK_SIZE)
    elseif nbytes + leafletbytes < BLOCK_SIZE
        GerminatingVine{rate}(
            ingest(trunk, leaflet), leaf, nbytes + leafletbytes)
    else
        trunkfill = (BLOCK_SIZE - nbytes) ÷ sizeof(U)
        trunk = convert(Sponge, ingest(trunk, view(leaflet, 1:trunkfill)))
        vine = Vine(ingest(trunk, K12_ZEROBLOCK_SUFFIX), leaf, UInt(BLOCK_SIZE))
        ingest(vine, view(leaflet, trunkfill+1:length(leaflet)))
    end
end

function ingest((; trunk, leaf, nbytes)::GerminatingVine{rate}, x::U) where {rate, U<:UInt8to64}
    if nbytes + sizeof(U) < BLOCK_SIZE
        GerminatingVine{rate}(ingest(trunk, x), leaf, nbytes + sizeof(U))
    elseif nbytes + sizeof(U) == BLOCK_SIZE
        Vine{rate}(convert(Sponge, ingest(trunk, x)),
                   leaf, nbytes + sizeof(U))
    else
        vine = GerminatingVine(trunk, leaf, nbytes)
        for byte in ntupleinterpret(UInt8, x)
            vine = ingest(vine, byte)
        end
        vine
    end
end

function ingest((; trunk, leaf, nbytes)::Vine{rate}, leaflet::AbstractVector{U}) where {rate, U<:UInt8to64}
    partial_bytes = (nbytes % BLOCK_SIZE) + length(leaflet)*sizeof(U)
    total_bytes = nbytes + length(leaflet)*sizeof(U)
    if partial_bytes < BLOCK_SIZE
        Vine(trunk, ingest(leaf, leaflet), total_bytes)
    elseif nbytes % sizeof(U) == 0
        ntofill = (BLOCK_SIZE - (nbytes % BLOCK_SIZE)) ÷ sizeof(U)
        leaf = ingest(leaf, view(leaflet, 1:ntofill))
        cv = squeeze(NTuple{4, UInt64}, pad(leaf, K12_SUFFIXES.leaf))
        trunk = ingest(trunk, cv)
        leaf = ByteSponge{rate}()
        for block in partition(ntofill+1:length(leaflet), BLOCK_SIZE÷sizeof(U))
            if length(block) == BLOCK_SIZE÷sizeof(U)
                trunk = ingest(trunk, NTuple{4, UInt64}, reinterpret(UInt64, view(leaflet, block)))
            else
                leaf = ingest(ByteSponge{rate}(), view(leaflet, block))
            end
        end
        Vine(trunk, leaf, total_bytes)
    else
        (@noinline () -> @error "Oh no! Something has gone very wrong, this message should never be printed 😟")()
        Vine(trunk, leaf, nbytes)
    end
end

function ingest(vine::Vine{rate}, x::U) where {rate, U<:UInt8to64}
    if (vine.nbytes % BLOCK_SIZE) + sizeof(U) < BLOCK_SIZE
        Vine(vine.trunk, ingest(vine.leaf, x), vine.nbytes + sizeof(U))
    elseif (vine.nbytes % BLOCK_SIZE) + sizeof(U) == BLOCK_SIZE
        Vine(ingest(vine.trunk, squeeze(NTuple{4, UInt64}, vine.leaf)),
             ingest(ByteSponge{rate}(), x), vine.nbytes + sizeof(U))
    else # `x` is not aligned with the partial data, so ingest byte by byte
        for byte in ntupleinterpret(UInt8, x)
            vine = ingest(vine, byte)
        end
        vine
    end
end

"""
    ingest_length(accumulator::Union{<:ByteSponge, <:AbstractVine}, val::UInt, ::Val{bufsize}=Val{8}())

Ingest a right-encoded form of `val` into `accumulator`, allowing the encoding to be up to `bufsize` bytes.
"""
function ingest_length(accum::Union{<:ByteSponge, <:AbstractVine}, val::UInt, ::Val{bufsize}=Val{8}()) where {bufsize}
    buffer = ntuple(_ -> 0x00, Val{bufsize}())
    point = 0
    while (val > 0)
        buffer = setindex(buffer, UInt8(val % 2^8), point+=1)
        val ÷= 2^8
    end
    for i in point:-1:1
        accum = ingest(accum, buffer[i])
    end
    ingest(accum, UInt8(point))
end

function ingest_length(accum::Union{<:ByteSponge, <:AbstractVine}, arr::AbstractVector{U}, ::Val{bufsize}=Val{8}()) where {bufsize, U <: UInt8to64}
    ingest_length(accum, UInt(length(arr)) * sizeof(U), Val{bufsize}())
end

"""
    finalise(vine::AbstractVine) -> ByteSponge

Finalise `vine` by folding in the current leaf returning the `pad`ded trunk.
"""
function finalise((; trunk, leaf, nbytes)::AbstractVine{rate}) where {rate}
    if leaf.byte != 1
        leaf = pad(leaf, K12_SUFFIXES.leaf)
        cv = squeeze(NTuple{4, UInt64}, leaf)
        trunk = ingest(trunk, cv)
    end
    sponge = convert(ByteSponge, trunk)
    if nbytes <= BLOCK_SIZE
        pad(sponge, K12_SUFFIXES.one)
    else
        sponge = ingest_length(sponge, nbytes ÷ BLOCK_SIZE)
        sponge = ingest(sponge, 0xffff)
        pad(sponge, K12_SUFFIXES.many)
    end
end
