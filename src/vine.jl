"""
    AbstractCoralVine{rate}

A binary tree where at each node, only one child itself has children.

```text
trunk ◆
     / \\
    ◆   ◇ leaf
   / \\
  ◆   ◇ leaf
  ┆
```

This esentially creates a "vine" like structure, with a single "trunk" that has
terminal "leaves" sprouting off it.

The trunk and its leaves are [sponges](@ref AbstractSponge), and it is this
"vine of sponges" conception that lends itself to the "Coral Vine" name (since
coral sponges are the only sponge that grows I'm aware of).

Since the final result of the vine is obtained from folding all leaves into it sequentially,
we can represent the intermediate state of the vine with just the "trunk" and
the most recent leaf (having folded all completed leaves into the trunk).

```text
                      ╭──────╮
                      │ Leaf │
╭─────────╮      □    ╰─┬────╯
│  Trunk  ├──┬───┴───┬──╯
╰─────────╯  □       □
```

However, at the very start data is absorbed into the trunk itself. After it is
full, all subsequent data goes to the leaves, and the trunk only absorbs values
squeezed from leaves.

The distinction between these two stages is made with the subtypes
[`CoralVineSeedling`](@ref) and [`CoralVine`](@ref).
"""
abstract type AbstractCoralVine{rate} end

"""
    CoralVineSeedling{rate}

A `vine` that hasn't started to grow leaves yet (see the [`AbstractCoralVine`](@ref) docstring).

```text
╭─────────╮
│  Trunk  ├───
╰─────────╯
```

See also: [`CoralVine`](@ref), [`absorb`](@ref), [`absorb_length`](@ref), and [`finalise`](@ref).
"""
struct CoralVineSeedling{rate} <: AbstractCoralVine{rate}
    trunk::ByteSponge{rate}
    nbytes::UInt
end

CoralVineSeedling{rate}() where {rate} =
    CoralVineSeedling{rate}(ByteSponge{rate}(), 0)

"""
    CoralVine{rate}

A `vine` with leaves (see the [`k12`](@ref) docstring).

```text
                      ╭──────╮
                      │ Leaf │
╭─────────╮      □    ╰─┬────╯
│  Trunk  ├──┬───┴───┬──╯
╰─────────╯  □       □
```

See also: [`CoralVineSeedling`](@ref), [`absorb`](@ref), [`absorb_length`](@ref), and [`finalise`](@ref).
"""
struct CoralVine{rate} <: AbstractCoralVine{rate}
    trunk::Sponge{rate}
    leaf::ByteSponge{rate}
    nbytes::UInt
end

"""
    absorb(vine::AbstractCoralVine, leaflet::AbstractVector{<:Unsigned})
    absorb(vine::AbstractCoralVine, x::Unsigned)

Ingest `leaflet`/`x` into `vine`, this may go into the leaves or the trunk depending
on the vine type, and may cause the current leaf to be folded into the trunk,
and a new leaf grown.
"""
function absorb((; trunk, nbytes)::CoralVineSeedling{rate}, leaflet::AbstractVector{U}) where {rate, U<:Unsigned}
    leafletbytes = length(leaflet) * sizeof(U)
    if nbytes == 0 && leafletbytes == BLOCK_SIZE
        zerostate = absorb(EMPTY_STATE, Val{rate2cap(rate)}(), reinterpret(UInt64, leaflet))
        CoralVine{rate}(absorb(Sponge{rate}(zerostate, 1 + (BLOCK_SIZE ÷ 8) % rate), K12_ZEROBLOCK_SUFFIX),
                        ByteSponge{rate}(), BLOCK_SIZE)
    elseif nbytes + leafletbytes < BLOCK_SIZE
        CoralVineSeedling{rate}(
            absorb(trunk, leaflet), nbytes + leafletbytes)
    else
        trunkfill = (BLOCK_SIZE - nbytes) ÷ sizeof(U)
        trunk = convert(Sponge, absorb(trunk, view(leaflet, 1:trunkfill)))
        vine = CoralVine(absorb(trunk, K12_ZEROBLOCK_SUFFIX), ByteSponge{rate}(), UInt(BLOCK_SIZE))
        absorb(vine, view(leaflet, trunkfill+1:length(leaflet)))
    end
end

function absorb((; trunk, nbytes)::CoralVineSeedling{rate}, x::U) where {rate, U<:UInt8to64}
    if nbytes + sizeof(U) < BLOCK_SIZE
        CoralVineSeedling{rate}(absorb(trunk, x), nbytes + sizeof(U))
    elseif nbytes + sizeof(U) == BLOCK_SIZE
        CoralVine{rate}(convert(Sponge, absorb(trunk, x)), ByteSponge{rate}(), nbytes + sizeof(U))
    else
        vine = CoralVineSeedling(trunk, nbytes)
        for byte in ntupleinterpret(UInt8, x)
            vine = absorb(vine, byte)
        end
        vine
    end
end

function absorb((; trunk, leaf, nbytes)::CoralVine{rate}, leaflet::AbstractVector{U}) where {rate, U<:UInt8to64}
    partial_bytes = (nbytes % BLOCK_SIZE) + length(leaflet)*sizeof(U)
    total_bytes = nbytes + length(leaflet)*sizeof(U)
    if partial_bytes < BLOCK_SIZE
        CoralVine(trunk, absorb(leaf, leaflet), total_bytes)
    elseif nbytes % sizeof(U) == 0
        ntofill = (BLOCK_SIZE - (nbytes % BLOCK_SIZE)) ÷ sizeof(U)
        leaf = absorb(leaf, view(leaflet, 1:ntofill))
        cv = squeeze(NTuple{4, UInt64}, pad(leaf, K12_SUFFIXES.leaf))
        trunk = absorb(trunk, cv)
        leaf = ByteSponge{rate}()
        for block in partition(ntofill+1:length(leaflet), BLOCK_SIZE÷sizeof(U))
            if length(block) == BLOCK_SIZE÷sizeof(U)
                trunk = absorb(trunk, NTuple{4, UInt64}, reinterpret(UInt64, view(leaflet, block)))
            else
                leaf = absorb(ByteSponge{rate}(), view(leaflet, block))
            end
        end
        CoralVine(trunk, leaf, total_bytes)
    else
        (@noinline () -> @error "Oh no! Something has gone very wrong, this message should never be printed 😟")()
        CoralVine(trunk, leaf, nbytes)
    end
end

function absorb(vine::CoralVine{rate}, x::U) where {rate, U<:UInt8to64}
    if (vine.nbytes % BLOCK_SIZE) + sizeof(U) < BLOCK_SIZE
        CoralVine(vine.trunk, absorb(vine.leaf, x), vine.nbytes + sizeof(U))
    elseif (vine.nbytes % BLOCK_SIZE) + sizeof(U) == BLOCK_SIZE
        CoralVine(absorb(vine.trunk, squeeze(NTuple{4, UInt64}, vine.leaf)),
             absorb(ByteSponge{rate}(), x), vine.nbytes + sizeof(U))
    else # `x` is not aligned with the partial data, so absorb byte by byte
        for byte in ntupleinterpret(UInt8, x)
            vine = absorb(vine, byte)
        end
        vine
    end
end

"""
    absorb_length(accumulator::Union{<:ByteSponge, <:AbstractCoralVine}, val::UInt, ::Val{bufsize}=Val{8}())

Ingest a right-encoded form of `val` into `accumulator`, allowing the encoding to be up to `bufsize` bytes.
"""
function absorb_length(accum::Union{<:ByteSponge, <:AbstractCoralVine}, val::UInt, ::Val{bufsize}=Val{8}()) where {bufsize}
    buffer = ntuple(_ -> 0x00, Val{bufsize}())
    point = 0
    while (val > 0)
        buffer = setindex(buffer, UInt8(val % 2^8), point+=1)
        val ÷= 2^8
    end
    for i in point:-1:1
        accum = absorb(accum, buffer[i])
    end
    absorb(accum, UInt8(point))
end

function absorb_length(accum::Union{<:ByteSponge, <:AbstractCoralVine}, arr::AbstractVector{U}, ::Val{bufsize}=Val{8}()) where {bufsize, U <: UInt8to64}
    absorb_length(accum, UInt(length(arr)) * sizeof(U), Val{bufsize}())
end

"""
    finalise(vine::AbstractCoralVine) -> ByteSponge

Finalise `vine` by folding in the current leaf returning the `pad`ded trunk.
"""
function finalise(vine::AbstractCoralVine{rate}) where {rate}
    (; trunk, nbytes) = vine
    if vine isa CoralVine && vine.leaf.byte != 0
        leaf = pad(vine.leaf, K12_SUFFIXES.leaf)
        cv = squeeze(NTuple{4, UInt64}, leaf)
        trunk = absorb(trunk, cv)
    end
    sponge = convert(ByteSponge, trunk)
    if nbytes <= BLOCK_SIZE
        pad(sponge, K12_SUFFIXES.one)
    else
        sponge = absorb_length(sponge, nbytes ÷ BLOCK_SIZE)
        sponge = absorb(sponge, 0xffff)
        pad(sponge, K12_SUFFIXES.many)
    end
end
