function k12_singlethreaded(message::AbstractVector{<:UInt8to64}, customisation::AbstractVector{<:UInt8to64})
    vine = absorb(CoralVineSeedling{RATE}(), message)
    vine = absorb(vine, customisation)
    vine = absorb_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function k12_singlethreaded(io::IO, customisation::AbstractVector{<:UInt8to64})
    block = Vector{UInt8}(undef, BLOCK_SIZE)
    vine = CoralVineSeedling{RATE}()
    while (size = readbytes!(io, block)) > 0
        vine = absorb(vine, view(block, 1:size))
    end
    vine = absorb(vine, customisation)
    vine = absorb_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function k12_singlethreaded_simd(message::AbstractVector{U}, customisation::AbstractVector{<:UInt8to64}) where {U<:UInt8to64}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return k12_singlethreaded(message, customisation)
    end
    slices = slice_message(U, length(message), Val{SIMD_FACTOR}())
    vine = absorb(CoralVineSeedling{RATE}(), view(message, slices.init))::CoralVine{RATE}
    bsize = BLOCK_SIZE÷sizeof(U)
    for simd_start in slices.blocks
        blocks = ntuple(i -> reinterpret(UInt64, view(message, simd_start+(i-1)*bsize:simd_start+i*bsize-1)),
                        Val{SIMD_FACTOR}())
        u64x4x4 = turboshake(NTuple{4, UInt64}, blocks, K12_SUFFIXES.leaf)
        trunk = absorb(vine.trunk, ntupleinterpret(UInt64, u64x4x4))
        vine = CoralVine(trunk, vine.leaf, vine.nbytes + BLOCK_SIZE * SIMD_FACTOR)
    end
    vine = absorb(vine, view(message, slices.tail))
    vine = absorb(vine, customisation)
    vine = absorb_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

function slice_message(U::Type, mlen::Int, ::Val{block_multiple}) where {block_multiple}
    u_block_size = BLOCK_SIZE÷sizeof(U)
    m_block_size = block_multiple * u_block_size
    init = if mlen >= u_block_size 1:min(mlen, u_block_size) else 1:0 end
    blocks_end = ((mlen - u_block_size) ÷ m_block_size) * m_block_size + last(init)
    blocks = last(init)+1:m_block_size:blocks_end-m_block_size+1
    tail = blocks_end+1:mlen
    (; init, blocks, tail)
end

function k12_multithreaded(message::AbstractVector{U}, customisation::AbstractVector{<:UInt8to64}) where {U<:UInt8to64}
    if length(message) <= BLOCK_SIZE÷sizeof(U)
        return k12_singlethreaded(message, customisation)
    end
    slices = slice_message(U, length(message), Val{1}())
    cvs = Vector{UInt64}(undef, 4 * length(slices.blocks))
    @threads for i in 1:length(slices.blocks)
        bstart = slices.blocks[i]
        block = reinterpret(UInt64, view(message, bstart:bstart-1+BLOCK_SIZE÷sizeof(U)))
        cv = @inline turboshake(NTuple{4, UInt64}, block, K12_SUFFIXES.leaf)
        for j in 1:4
            cvs[4(i-1)+j] = cv[j]
        end
    end
    vine = absorb(CoralVineSeedling{RATE}(), view(message, slices.init))::CoralVine{RATE}
    vine = CoralVine(absorb(vine.trunk, cvs), vine.leaf, vine.nbytes + length(cvs) * BLOCK_SIZE ÷ 4)
    vine = absorb(vine, view(message, slices.tail))
    vine = absorb(vine, customisation)
    vine = absorb_length(vine, customisation)
    squeeze(UInt128, finalise(vine))
end

"""
    k12(data::Union{IO, String, AbstractVector{<:Unsigned}},
        customisation::Union{String, AbstractVector{<:Unsigned}};
        thread::Bool=true) -> UInt128

Hash `data` with the KangarooTwelve scheme and a `customisation` value, using
multithreading when `thread` is `true`.

This scheme presents a good balance of *Simplicity*, *Security*, and *Speed*.

# Extended help

The KangarooTwelve hashing scheme works by splitting the input data into ``n``
$BLOCK_SIZE-byte blocks (``S₀``, ``S₁``, …, ``Sₙ₋₁``) which are individually
hashed with [`TurboSHAKE128`](@ref turboshake) to produce 32-byte "chaining
values" (CVs), which are put together and absorbed to produce the final state.

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
function k12(data::AbstractVector{<:Unsigned}, customisation::AbstractVector{<:Unsigned}=UInt8[]; thread::Bool=true)
    if thread
        k12_multithreaded(data, customisation)
    else
        k12_singlethreaded(data, customisation)
    end
end

function k12(data::IO, customisation::AbstractVector{<:Unsigned}=UInt8[]; thread::Bool=true)
    thread || return k12_singlethreaded(data, customisation)
    try
        buf = Mmap.mmap(data)
        k12(buf, customisation)
    catch _
        k12_singlethreaded(data, customisation)
    end
end

k12(s::AbstractString, customisation::AbstractVector{<:Unsigned}=UInt8[]; thread::Bool=true) =
    k12(codeunits(String(s)), customisation; thread)

k12(data::Union{<:AbstractVector{<:Unsigned}, <:IO, <:AbstractString}, customisation::AbstractString) =
    k12(data, codeunits(String(customisation)))
