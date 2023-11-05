function k12_singlethreaded(message::AbstractVector{<:UInt8to64}, customisation::AbstractVector{<:UInt8to64})
    vine = ingest(GerminatingVine{RATE}(), message)
    vine = ingest(vine, customisation)
    vine = ingest_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function k12_singlethreaded(io::IO, customisation::AbstractVector{<:UInt8to64})
    block = Vector{UInt8}(undef, BLOCK_SIZE)
    vine = GerminatingVine{RATE}()
    while (size = readbytes!(io, block)) > 0
        vine = ingest(vine, view(block, 1:size))
    end
    vine = ingest(vine, customisation)
    vine = ingest_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function k12_singlethreaded_simd(message::AbstractVector{U}, customisation::AbstractVector{<:UInt8to64}) where {U<:UInt8to64}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return k12_singlethreaded(message, customisation)
    end
    slices = slice_message(U, length(message), Val{SIMD_FACTOR}())
    vine = ingest(GerminatingVine{RATE}(), view(message, slices.init))::Vine{RATE}
    bsize = BLOCK_SIZE÷sizeof(U)
    for simd_start in slices.blocks
        blocks = ntuple(i -> reinterpret(UInt64, view(message, simd_start+(i-1)*bsize:simd_start+i*bsize-1)),
                        Val{SIMD_FACTOR}())
        u64x4x4 = turboshake(NTuple{4, UInt64}, blocks, K12_SUFFIXES.leaf)
        trunk = ingest(vine.trunk, ntupleinterpret(UInt64, u64x4x4))
        vine = Vine(trunk, vine.leaf, vine.nbytes + BLOCK_SIZE * SIMD_FACTOR)
    end
    vine = ingest(vine, view(message, slices.tail))
    vine = ingest(vine, customisation)
    vine = ingest_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function slice_message(U::Type, mlen::Int, ::Val{block_multiple}) where {block_multiple}
    u_block_size = BLOCK_SIZE÷sizeof(U)
    m_block_size = block_multiple * u_block_size
    init = if mlen >= u_block_size 1:min(mlen, u_block_size) else 1:0 end
    blocks_end = if block_multiple > 0
        ((mlen - u_block_size) ÷ m_block_size) * m_block_size + last(init)
    else last(init) end
    blocks = last(init)+1:m_block_size:blocks_end-m_block_size+1
    tail = blocks_end+1:mlen
    (; init, blocks, tail)
end

function k12_multithreaded(message::AbstractVector{U}, customisation::AbstractVector{<:UInt8to64}) where {U<:UInt8to64}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return k12_singlethreaded(message, customisation)
    end
    slices = slice_message(U, length(message), Val{1}())
    blocks = [reinterpret(UInt64, view(message, bstart:bstart-1+BLOCK_SIZE÷sizeof(U)))
              for bstart in slices.blocks]
    cvs = Vector{UInt64}(undef, 4 * length(blocks))
    @threads for i in 1:length(blocks)
        block = blocks[i]
        cv = turboshake(NTuple{4, UInt64}, block, K12_SUFFIXES.leaf)
        for j in 1:4
            cvs[4(i-1)+j] = cv[j]
        end
    end
    vine = ingest(GerminatingVine{RATE}(), view(message, slices.init))
    for cv in cvs
        trunk = ingest(vine.trunk, cv)
        vine = Vine(trunk, vine.leaf, vine.nbytes + BLOCK_SIZE ÷ 4)
    end
    vine = ingest(vine, view(message, slices.tail))
    vine = ingest(vine, customisation)
    vine = ingest_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

"""
    k12(data::Union{IO, String, AbstractVector{Unsigned}}) -> UInt128

Hash `data` with the KangarooTwelve scheme.

This scheme presents a good balance of *Simplicity*, *Security*, and *Speed*.

# Extended help

The KangarooTwelve hashing scheme works by splitting the input data into ``n``
$BLOCK_SIZE-byte blocks (``S₀``, ``S₁``, …, ``Sₙ₋₁``) which are individually
hashed with [`TurboSHAKE128`](@ref turboshake) to produce 32-byte "chaining
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
SIMD and multithreading.
"""
function k12(data::AbstractVector{<:Unsigned}, customisation::AbstractVector{<:Unsigned}=UInt8[])
    if true # length(data) < heuristic
        k12_singlethreaded(data, customisation)
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
