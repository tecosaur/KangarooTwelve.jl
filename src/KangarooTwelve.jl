module KangarooTwelve

using Base.Threads
using SIMD
using Mmap

# See <https://www.fress.io/story/k12>
# For SIMD oppotunities, see <https://docs.rs/keccak/latest/src/keccak/lib.rs.html>.
# With native code, <https://lists.gnupg.org/pipermail/gcrypt-devel/2022-July/005360.html>
# and <https://github.com/XKCP/XKCP/blob/master/lib/low/KeccakP-1600/AVX2/KeccakP-1600-AVX2.s>
# for some hints as to what good assembly should look like.
#
# The performance should be similar to <https://keccak.team/sw_performance.html>

export k12

import Base.Cartesian.@ntuple
import Base.Iterators.partition
import Base.setindex

const UInt8to64 = Union{UInt64, UInt32, UInt16, UInt8}

const CAPACITY = 256
const RATE = 25 - CAPACITY ÷ 64
const BLOCK_SIZE = 8192

const SIMD_FACTOR = 4

const K12_SUFFIXES = (one=0x07, many=0x06, leaf=0x0b)
const K12_ZEROBLOCK_SUFFIX = 0x0000000000000003

include("keccak.jl")
include("turboshake.jl")
include("sponge.jl")
include("vine.jl")
include("k12.jl")

include("throughput.jl")

end
