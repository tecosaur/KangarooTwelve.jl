function k12_singlethreaded(message::AbstractVector{<:UInt8to64}, customisation::AbstractVector{<:UInt8to64})
    vine = ingest(GerminatingVine{RATE}(), message)
    vine = ingest(vine, customisation)
    vine = ingest_length(vine, customisation)
    vine = finalise(vine)
    squeeze(UInt128, vine)
end

const SIMD_BLOCKS = ntuple(i -> 1+(i-1)*BLOCK_SIZE÷sizeof(UInt64):i*BLOCK_SIZE÷sizeof(UInt64),
                           Val{SIMD_FACTOR}())

function k12_singlethreaded_simd(message::AbstractVector{U}, customisation::AbstractVector{V}) where {U <: UInt8to64, V <: UInt8to64}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return k12_singlethreaded(message, customisation)
    end
    slices = slice_message(message, customisation, Val{SIMD_FACTOR}())
    n, sponge = 0, Sponge(slices.init)
    for sblock in slices.simd
        blocks = ntuple(i -> view(sblock, SIMD_BLOCKS[i]), Val{SIMD_FACTOR}())
        out_4u32 = turboshake(UInt32, blocks, K12_SUFFIXES.leaf)
        out1u64, out2u64 = reinterpret(NTuple{2, UInt64}, out_4u32)
        sponge, n = ingest(sponge, out1u64, out2u64), n + SIMD_FACTOR
    end
    for block in slices.blocks
        sponge, n = ingest(sponge, UInt32, block), n+1
    end
    if !isempty(slices.tail)
        sponge, n = ingest(sponge, slices.tail), n+1
    end
    sponge = ingest(ingest_length(sponge, n), 0xffff, 0x01)
    sponge = pad(sponge, K12_SUFFIXES.many)
    squeeze(UInt128, sponge)
end

function slice_message((U, mlen)::Tuple{Type, Int}, (V, clen)::Tuple{Type, Int}, ::Val{simd_factor}) where {simd_factor}
    u_block_size = BLOCK_SIZE÷sizeof(U)
    v_block_size = BLOCK_SIZE÷sizeof(V)
    simd_block_size = simd_factor * u_block_size
    init = if mlen >= u_block_size 1:min(mlen, u_block_size) else 1:0 end
    message_simd_end = if simd_factor > 0
        ((mlen - u_block_size) ÷ simd_block_size) * simd_block_size + last(init)
    else last(init) end
    tail = last(simd_region)+1:mlen
    (; init, simd_region, tail)
end

function slice_message(message::AbstractVector{U}, customisation::AbstractVector{V}, ::Val{simd_factor}) where {U, V, simd_factor}
    slices = slice_message((U, length(message)), (V, length(customisation)), Val{simd_factor}())
    (; init = view(message, slices.init),
     simd = if simd_factor > 0
         partition(uinterpret(UInt64, view(message, slices.simd_region)), BLOCK_SIZE÷sizeof(UInt64) * simd_factor)
     else () end,
     blocks = partition(uinterpret(UInt64, view(message, slices.block_region)), BLOCK_SIZE÷sizeof(UInt64)),
     tail = view(message, slices.tail),
     custom_neck = view(customisation, slices.custom_neck),
     custom_blocks = partition(uinterpret(UInt64, view(customisation, slices.custom_block_region)), BLOCK_SIZE÷sizeof(UInt64)),
     custom_tail = view(customisation, slices.custom_tail))
end

# TODO: Better multithreaded implementation + heuristic
function k12_multithreaded(message::AbstractVector{U}) where {U <: UInt8to64}
    chunks = partition(message, SIMD_FACTOR * BLOCK_SIZE÷sizeof(U))
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
    n, sponge = 0, Sponge(block)
    while (size = readbytes!(io, block)) > 0
        out32 = turboshake(UInt32, view(block, 1:size), K12_SUFFIXES.leaf)
        sponge, n = ingest(sponge, out32), n+1
    end
    sponge = ingest(ingest_length(sponge, n), 0xffff, 0x01)
    sponge = pad(sponge, K12_SUFFIXES.many)
    squeeze(UInt128, sponge)
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
