var documenterSearchIndex = {"docs":
[{"location":"benchmark/#Benchmark","page":"Benchmark","title":"Benchmark","text":"","category":"section"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"The single-threaded and multi-threaded variants of KangarooTwelve from this package are compared against:","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"crc32c from the CRC32c standard library\nsha256sum from the SHA standard library\nmd5 from MD5.jl\nblake3 from Blake3Hash.jl","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"(Image: image)","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"System information","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"Intel i5-13600k\nLinux 6.6 (openSUSE)","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"Version Information","category":"page"},{"location":"benchmark/","page":"Benchmark","title":"Benchmark","text":"Julia 1.10-rc1\nSHA v0.7.0\nMD5 v0.2.2\nBlake3Hash v0.3.0","category":"page"},{"location":"#KangarooTwelve","page":"KangarooTwelve","title":"KangarooTwelve","text":"","category":"section"},{"location":"#Background","page":"KangarooTwelve","title":"Background","text":"","category":"section"},{"location":"","page":"KangarooTwelve","title":"KangarooTwelve","text":"A pure-Julia implementation of the KangarooTwelve hashing scheme, so named because it consists of 12 rounds of Keccak (TurboSHAKE128) with kangaroo hopping, allowing for parallel tree-ish hashing termed \"leaves stapled to a pole\".","category":"page"},{"location":"","page":"KangarooTwelve","title":"KangarooTwelve","text":"This scheme presents a particularly good balance of:","category":"page"},{"location":"","page":"KangarooTwelve","title":"KangarooTwelve","text":"Simplicity (Keccak + sponge + hopping)\nSecurity (128-bit)\nSpeed (up to ~2bytes/cycle)","category":"page"},{"location":"","page":"KangarooTwelve","title":"KangarooTwelve","text":"It is currently an IETF draft.","category":"page"},{"location":"#Usage","page":"KangarooTwelve","title":"Usage","text":"","category":"section"},{"location":"","page":"KangarooTwelve","title":"KangarooTwelve","text":"k12","category":"page"},{"location":"#KangarooTwelve.k12","page":"KangarooTwelve","title":"KangarooTwelve.k12","text":"k12(data::Union{IO, String, AbstractVector{<:Unsigned}},\n    customisation::Union{String, AbstractVector{<:Unsigned}};\n    thread::Bool=true) -> UInt128\n\nHash data with the KangarooTwelve scheme and a customisation value, using multithreading when thread is true.\n\nThis scheme presents a good balance of Simplicity, Security, and Speed.\n\nExtended help\n\nThe KangarooTwelve hashing scheme works by splitting the input data into n 8192-byte blocks (S₀, S₁, …, Sₙ₁) which are individually hashed with TurboSHAKE128 to produce 32-byte \"chaining values\" (CVs), which are put together and absorbed to produce the final state.\n\n               ╭────╮ ╭────╮   ╭────╮ ╭────╮\n               │ S₁ │ │ S₂ │   │Sₙ₋₂│ │Sₙ₋₁│\n               ╰─┬──╯ ╰─┬──╯   ╰─┬──╯ ╰─┬──╯\n                 │110   │110     │110   │110\n                 ▼      ▼        ▼      ▼\n╭────────╮110⁶²╭─┴──╮ ╭─┴──╮   ╭─┴──╮ ╭─┴──╮\n│   S₀   ├─────┤ CV ├─┤ CV ├╴╍╶┤ CV ├─┤ CV ├─(n-1)(FFFF)(01)──▶─┤HASH│\n╰────────╯     ╰────╯ ╰────╯   ╰────╯ ╰────╯\n\nThis scheme has been described as \"leaves on a vine\". The hashing of blocks S₁ to Sₙ₁ is embarassingly parallel, and can be accelerated with both SIMD and multithreading.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Internals","page":"Internals","title":"Internals","text":"","category":"section"},{"location":"internals/#Keccak","page":"Internals","title":"Keccak","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"KangarooTwelve.keccak_p1600","category":"page"},{"location":"internals/#KangarooTwelve.keccak_p1600","page":"Internals","title":"KangarooTwelve.keccak_p1600","text":"keccak_p1600(state::NTuple{25, UInt64}, ::Val{nrounds}=Val{12}())\n\nApply the Keccak-p[nrounds, 1600] permutation to state. This is formally defined in the Keccak reference  and formalised in FIPS 202.\n\nExtended help\n\nThis is a variable-round permutation, with up to 24 rounds, where the last round performed matches the final Keccak-f permutation round.\n\nEach round of permutation consists of five steps, termed θ, ρ, π, χ, and ι. These all operate on a 200-byte state.\n\nThese steps view the state in a number of different configurations.\n\n             ┌─┬─┬─┬─┬─┐\n            ┌─┬─┬─┬─┬─┐┤\n           ┌─┬─┬─┬─┬─┐┤┤\n          ┌─┬─┬─┬─┬─┐┤┤┤\n         ┌─┬─┬─┬─┬─┐┤┤┤┤\n        ┌─┬─┬─┬─┬─┐┤┤┤┤┘\n       ┌─┬─┬─┬─┬─┐┤┤┤┤┘\n      ┌─┬─┬─┬─┬─┐┤┤┤┤┘\n      ├─┼─┼─┼─┼─┤┤┤┤┘\n      ├─┼─┼─┼─┼─┤┤┤┘\n      ├─┼─┼─┼─┼─┤┤┘        ┌─┐ bit\n      ├─┼─┼─┼─┼─┤┘         └─┘         ┌─┐\n      └─┴─┴─┴─┴─┘                     ┌─┐┘\n         state                       ┌─┐┘\n                                    ┌─┐┘\n                   ┌─┐             ┌─┐┘\n                   ├─┤ column     ┌─┐┘\n     row           ├─┤           ┌─┐┘\n ┌─┬─┬─┬─┬─┐       ├─┤          ┌─┐┘  lane\n └─┴─┴─┴─┴─┘       ├─┤          └─┘\n                   └─┘\n\nθ step Compute the parity of five columns, and xor-diffuse their parity into nearby columns.\nρ step Bitwise-rotate each of the 25 lanes by a different triangular number.\nπ step Permute each lane in a fixed pattern.\nχ step Intra-row bitwise combination. This provides the non-linearity.\nι step The first lane is (xor-)mixed with a LFSR sequence across rounds. This serves to disrupt the symmetry of the scheme.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Turboshake","page":"Internals","title":"Turboshake","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"KangarooTwelve.turboshake","category":"page"},{"location":"internals/#KangarooTwelve.turboshake","page":"Internals","title":"KangarooTwelve.turboshake","text":"turboshake(output::Type, message::AbstractVector{<:UInt8to64},\n           delimsuffix::UInt8=0x80, ::Val{capacity} = Val{256}())\n\nProduce an output (Unsigned or NTuple{n, <:Unsigned}) value, by performing TurboSHAKE on message with a certain capacity and delimsuffix.\n\nSee also: absorb, pad, and squeeze.\n\n\n\n\n\n","category":"function"},{"location":"internals/#Sponge-and-Vine","page":"Internals","title":"Sponge & Vine","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"KangarooTwelve.AbstractSponge\nKangarooTwelve.Sponge\nKangarooTwelve.ByteSponge\nKangarooTwelve.AbstractCoralVine\nKangarooTwelve.CoralVineSeedling\nKangarooTwelve.CoralVine\nKangarooTwelve.absorb_length\nKangarooTwelve.finalise","category":"page"},{"location":"internals/#KangarooTwelve.AbstractSponge","page":"Internals","title":"KangarooTwelve.AbstractSponge","text":"AbstractSponge\n\nA sponge is a generalisation of hash functions and stream ciphers. It is in essance an iterated construction of a variable-input variable-output function based on a fixed length permutation and an internal state.\n\nThe way the sponge operates can be illustrated well with a diagram following the state.\n\n Message (padded)\n  └───┬──────┬──────┬──────┐      ┊ ┌──────┬──────┬───▶ Output\n      │      │      │      │      ┊ │      │      │\n ┌─┐  │  ╭─╮ │  ╭─╮ │  ╭─╮ │  ╭─╮ ┆ │  ╭─╮ │  ╭─╮ │\n │ │  ▼  │ │ ▼  │ │ ▼  │ │ ▼  │ │ ┆ │  │ │ │  │ │ │\nr│0├──⨁─▶│ ├─⨁─▶│ ├─⨁─▶│ ├─⨁─▶│ ├─┼─┴─▶│ ├─┴─▶│ ├─┴─▶ ╍\n │ │     │f│    │f│    │f│    │f│ ┆    │f│    │f│\n ├─┤     │ │    │ │    │ │    │ │ ┆    │ │    │ │\nc│0├────▶│ ├───▶│ ├───▶│ ├───▶│ ├─┼───▶│ ├───▶│ ├───▶ ╍\n └─┘     ╰─╯    ╰─╯    ╰─╯    ╰─╯ ┆    ╰─╯    ╰─╯\n                        Absorbing ┆ Squeezing\n\nFirst, the input is padded with a reversible padding rule and the state is initialised to zero.\n\nIn the absorbing phase, the input is processed in blocks of r bits, and xor'd with the first r bits of the state (leaving the remaining c bits unchanged), with absorptions interlieved with applications of the function f.\nIn the squeezing phase, the first r bits of the state are returned as the output. Shoult more output be wanted, f can simply be applied to the state again.\n\nWhile this construction is defined for any fixed-length permutation f, however we specialise on keccak_p1600.\n\nIt is possible to either incrementally either one lane at a time, or one byte at a time. To handle these two cases we have the Sponge and ByteSponge subtypes.\n\n\n\n\n\n","category":"type"},{"location":"internals/#KangarooTwelve.Sponge","page":"Internals","title":"KangarooTwelve.Sponge","text":"Sponge{rate}\n\nA Keccak state that keeps track of the last lane updated.\n\nFor more information on the sponge construction see AbstractSponge.\n\nSee also: absorb, pad, and squeeze.\n\n\n\n\n\n","category":"type"},{"location":"internals/#KangarooTwelve.ByteSponge","page":"Internals","title":"KangarooTwelve.ByteSponge","text":"ByteSponge{rate}\n\nA Keccak state that keeps track of the last byte updated.\n\nFor more information on the sponge construction see AbstractSponge.\n\nSee also: absorb, pad, and squeeze.\n\n\n\n\n\n","category":"type"},{"location":"internals/#KangarooTwelve.AbstractCoralVine","page":"Internals","title":"KangarooTwelve.AbstractCoralVine","text":"AbstractCoralVine{rate}\n\nA binary tree where at each node, only one child itself has children.\n\ntrunk ◆\n     / \\\n    ◆   ◇ leaf\n   / \\\n  ◆   ◇ leaf\n  ┆\n\nThis esentially creates a \"vine\" like structure, with a single \"trunk\" that has terminal \"leaves\" sprouting off it.\n\nThe trunk and its leaves are sponges, and it is this \"vine of sponges\" conception that lends itself to the \"Coral Vine\" name (since coral sponges are the only sponge that grows I'm aware of).\n\nSince the final result of the vine is obtained from folding all leaves into it sequentially, we can represent the intermediate state of the vine with just the \"trunk\" and the most recent leaf (having folded all completed leaves into the trunk).\n\n                      ╭──────╮\n                      │ Leaf │\n╭─────────╮      □    ╰─┬────╯\n│  Trunk  ├──┬───┴───┬──╯\n╰─────────╯  □       □\n\nHowever, at the very start data is absorbed into the trunk itself. After it is full, all subsequent data goes to the leaves, and the trunk only absorbs values squeezed from leaves.\n\nThe distinction between these two stages is made with the subtypes CoralVineSeedling and CoralVine.\n\n\n\n\n\n","category":"type"},{"location":"internals/#KangarooTwelve.CoralVineSeedling","page":"Internals","title":"KangarooTwelve.CoralVineSeedling","text":"CoralVineSeedling{rate}\n\nA vine that hasn't started to grow leaves yet (see the AbstractCoralVine docstring).\n\n╭─────────╮\n│  Trunk  ├───\n╰─────────╯\n\nSee also: CoralVine, absorb, absorb_length, and finalise.\n\n\n\n\n\n","category":"type"},{"location":"internals/#KangarooTwelve.CoralVine","page":"Internals","title":"KangarooTwelve.CoralVine","text":"CoralVine{rate}\n\nA vine with leaves (see the k12 docstring).\n\n                      ╭──────╮\n                      │ Leaf │\n╭─────────╮      □    ╰─┬────╯\n│  Trunk  ├──┬───┴───┬──╯\n╰─────────╯  □       □\n\nSee also: CoralVineSeedling, absorb, absorb_length, and finalise.\n\n\n\n\n\n","category":"type"},{"location":"internals/#KangarooTwelve.absorb_length","page":"Internals","title":"KangarooTwelve.absorb_length","text":"absorb_length(accumulator::Union{<:ByteSponge, <:AbstractCoralVine}, val::UInt, ::Val{bufsize}=Val{8}())\n\nIngest a right-encoded form of val into accumulator, allowing the encoding to be up to bufsize bytes.\n\n\n\n\n\n","category":"function"},{"location":"internals/#KangarooTwelve.finalise","page":"Internals","title":"KangarooTwelve.finalise","text":"finalise(vine::AbstractCoralVine) -> ByteSponge\n\nFinalise vine by folding in the current leaf returning the padded trunk.\n\n\n\n\n\n","category":"function"},{"location":"internals/#State-operations","page":"Internals","title":"State operations","text":"","category":"section"},{"location":"internals/","page":"Internals","title":"Internals","text":"The absorb, pad, squeeze, and squeeze! are implemented for state tuples, sponges, and vines.","category":"page"},{"location":"internals/","page":"Internals","title":"Internals","text":"KangarooTwelve.absorb\nKangarooTwelve.pad\nKangarooTwelve.squeeze\nKangarooTwelve.squeeze!","category":"page"},{"location":"internals/#KangarooTwelve.absorb","page":"Internals","title":"KangarooTwelve.absorb","text":"absorb(state::NTuple{25, UInt64}, block::NTuple{rate, UInt64})\n\nIngest a single block of input (with implied rate) into state.\n\nThe first rate elements of the state are xor'd block, and then the state is permuted with keccak_p1600.\n\n\n\n\n\nabsorb(state::NTuple{25, UInt64}, ::Val{capacity}, message::AbstractVector{<:Unsigned})\n\nIngest message into state, with the rate calculated based on capacity.\n\nThis breaks message into rate-sized blocks and then absorbs them (as per absorb(state, ::NTuple{rate, UInt64})) in turn.\n\n\n\n\n\nabsorb(sponge::AbstractSponge, as::Type{<:Unsigned}, leaf::AbstractVector{<:Unsigned})\n\nIngest leaf into sponge by transforming it to an as via turboshake, and absorbing that result.\n\n\n\n\n\nabsorb(sponge::AbstractSponge, block::AbstractVector{<:Unsigned})\n\nIngest each element of block into sponge.\n\n\n\n\n\nabsorb(sponge::ByteSponge, x::Union{<:Unsigned, NTuple{N, <:Unsigned}})\n\nIngest x into sponge.\n\n\n\n\n\nabsorb(vine::AbstractCoralVine, leaflet::AbstractVector{<:Unsigned})\nabsorb(vine::AbstractCoralVine, x::Unsigned)\n\nIngest leaflet/x into vine, this may go into the leaves or the trunk depending on the vine type, and may cause the current leaf to be folded into the trunk, and a new leaf grown.\n\n\n\n\n\n","category":"function"},{"location":"internals/#KangarooTwelve.pad","page":"Internals","title":"KangarooTwelve.pad","text":"pad(state::NTuple{25, UInt64}, ::Val{capacity}, lastbyte::UInt, delimsuffix::UInt8)\n\nPerform \"padding\" of the lastbyte byte of state with delimsuffix.\n\nThis is the final modification to state in turboshake.\n\n\n\n\n\n","category":"function"},{"location":"internals/#KangarooTwelve.squeeze","page":"Internals","title":"KangarooTwelve.squeeze","text":"squeeze(outtype::Type, state::NTuple{25, UInt64}, ::Val{capacity})\n\nSqueeze an outtype out of state. outtype can be an Unsigned type or an unsigned NTuple.\n\n\n\n\n\nsqueeze(T::Type, sponge::AbstractSponge)\n\nSqueeze a T out of sponge.\n\n\n\n\n\n","category":"function"},{"location":"internals/#KangarooTwelve.squeeze!","page":"Internals","title":"KangarooTwelve.squeeze!","text":"squeeze!(output::Vector{UInt64}, state::NTuple{25, UInt64}, ::Val{capacity})\n\nSqueeze state into output.\n\n\n\n\n\nsqueeze!(output::Vector{UInt64}, sponge::AbstractSponge)\n\nSqueeze sponge into output.\n\n\n\n\n\n","category":"function"}]
}
