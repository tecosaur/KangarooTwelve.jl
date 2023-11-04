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
    for simd_start in slices.simd
        blocks = ntuple(i -> reinterpret(UInt64, view(message, simd_start+(i-1)*bsize:simd_start+i*bsize-1)),
                        Val{SIMD_FACTOR}())
        u64x4x4 = turboshake(NTuple{4, UInt64}, blocks, K12_SUFFIXES.leaf)
        trunk = ingest(vine.trunk, reinterpret(NTuple{16, UInt64}, u64x4x4))
        vine = Vine(trunk, vine.leaf, vine.nbytes + BLOCK_SIZE * SIMD_FACTOR)
    end
    vine = ingest(vine, view(message, slices.tail))
    vine = ingest(vine, customisation)
    vine = ingest_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function slice_message(U::Type, mlen::Int, ::Val{simd_factor}) where {simd_factor}
    u_block_size = BLOCK_SIZE÷sizeof(U)
    simd_block_size = simd_factor * u_block_size
    init = if mlen >= u_block_size 1:min(mlen, u_block_size) else 1:0 end
    simd_end = if simd_factor > 0
        ((mlen - u_block_size) ÷ simd_block_size) * simd_block_size + last(init)
    else last(init) end
    simd = last(init)+1:simd_block_size:simd_end-simd_block_size+1
    tail = simd_end+1:mlen
    (; init, simd, tail)
end

# TODO: Better multithreaded implementation + heuristic
function k12_multithreaded(message::AbstractVector{U}) where {U <: UInt8to64}
    chunks = partition(message, SIMD_FACTOR * BLOCK_SIZE÷sizeof(U))
    tasks = map(chunks) do chunk
        @spawn turboshake(UInt32, chunk, K12_SUFFIXES.leaf)
    end
    map(fetch, tasks)
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
