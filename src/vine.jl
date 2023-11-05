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

See also: [`CoralVine`](@ref), `ingest`, `ingest_length`, `finalise`.
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

See also: [`CoralVineSeedling`](@ref), `ingest`, `ingest_length`, `finalise`.
"""
struct CoralVine{rate} <: AbstractCoralVine{rate}
    trunk::Sponge{rate}
    leaf::ByteSponge{rate}
    nbytes::UInt
end

"""
    ingest(vine::AbstractCoralVine, leaflet::AbstractVector{<:Unsigned})
    ingest(vine::AbstractCoralVine, x::Unsigned)

Ingest `leaflet`/`x` into `vine`, this may go into the leaves or the trunk depending
on the vine type, and may cause the current leaf to be folded into the trunk,
and a new leaf grown.
"""
function ingest((; trunk, nbytes)::CoralVineSeedling{rate}, leaflet::AbstractVector{U}) where {rate, U<:Unsigned}
    leafletbytes = length(leaflet) * sizeof(U)
    if nbytes == 0 && leafletbytes == BLOCK_SIZE
        zerostate = ingest(EMPTY_STATE, Val{rate2cap(rate)}(), reinterpret(UInt64, leaflet))
        CoralVine{rate}(ingest(Sponge{rate}(zerostate, 1 + (BLOCK_SIZE ÷ 8) % rate), K12_ZEROBLOCK_SUFFIX),
                        ByteSponge{rate}(), BLOCK_SIZE)
    elseif nbytes + leafletbytes < BLOCK_SIZE
        CoralVineSeedling{rate}(
            ingest(trunk, leaflet), nbytes + leafletbytes)
    else
        trunkfill = (BLOCK_SIZE - nbytes) ÷ sizeof(U)
        trunk = convert(Sponge, ingest(trunk, view(leaflet, 1:trunkfill)))
        vine = CoralVine(ingest(trunk, K12_ZEROBLOCK_SUFFIX), ByteSponge{rate}(), UInt(BLOCK_SIZE))
        ingest(vine, view(leaflet, trunkfill+1:length(leaflet)))
    end
end

function ingest((; trunk, nbytes)::CoralVineSeedling{rate}, x::U) where {rate, U<:UInt8to64}
    if nbytes + sizeof(U) < BLOCK_SIZE
        CoralVineSeedling{rate}(ingest(trunk, x), nbytes + sizeof(U))
    elseif nbytes + sizeof(U) == BLOCK_SIZE
        CoralVine{rate}(convert(Sponge, ingest(trunk, x)), ByteSponge{rate}(), nbytes + sizeof(U))
    else
        vine = CoralVineSeedling(trunk, nbytes)
        for byte in ntupleinterpret(UInt8, x)
            vine = ingest(vine, byte)
        end
        vine
    end
end

function ingest((; trunk, leaf, nbytes)::CoralVine{rate}, leaflet::AbstractVector{U}) where {rate, U<:UInt8to64}
    partial_bytes = (nbytes % BLOCK_SIZE) + length(leaflet)*sizeof(U)
    total_bytes = nbytes + length(leaflet)*sizeof(U)
    if partial_bytes < BLOCK_SIZE
        CoralVine(trunk, ingest(leaf, leaflet), total_bytes)
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
        CoralVine(trunk, leaf, total_bytes)
    else
        (@noinline () -> @error "Oh no! Something has gone very wrong, this message should never be printed 😟")()
        CoralVine(trunk, leaf, nbytes)
    end
end

function ingest(vine::CoralVine{rate}, x::U) where {rate, U<:UInt8to64}
    if (vine.nbytes % BLOCK_SIZE) + sizeof(U) < BLOCK_SIZE
        CoralVine(vine.trunk, ingest(vine.leaf, x), vine.nbytes + sizeof(U))
    elseif (vine.nbytes % BLOCK_SIZE) + sizeof(U) == BLOCK_SIZE
        CoralVine(ingest(vine.trunk, squeeze(NTuple{4, UInt64}, vine.leaf)),
             ingest(ByteSponge{rate}(), x), vine.nbytes + sizeof(U))
    else # `x` is not aligned with the partial data, so ingest byte by byte
        for byte in ntupleinterpret(UInt8, x)
            vine = ingest(vine, byte)
        end
        vine
    end
end

"""
    ingest_length(accumulator::Union{<:ByteSponge, <:AbstractCoralVine}, val::UInt, ::Val{bufsize}=Val{8}())

Ingest a right-encoded form of `val` into `accumulator`, allowing the encoding to be up to `bufsize` bytes.
"""
function ingest_length(accum::Union{<:ByteSponge, <:AbstractCoralVine}, val::UInt, ::Val{bufsize}=Val{8}()) where {bufsize}
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

function ingest_length(accum::Union{<:ByteSponge, <:AbstractCoralVine}, arr::AbstractVector{U}, ::Val{bufsize}=Val{8}()) where {bufsize, U <: UInt8to64}
    ingest_length(accum, UInt(length(arr)) * sizeof(U), Val{bufsize}())
end

"""
    finalise(vine::AbstractCoralVine) -> ByteSponge

Finalise `vine` by folding in the current leaf returning the `pad`ded trunk.
"""
function finalise(vine::AbstractCoralVine{rate}) where {rate}
    (; trunk, nbytes) = vine
    if vine isa CoralVine && vine.leaf.byte != 1
        leaf = pad(vine.leaf, K12_SUFFIXES.leaf)
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
