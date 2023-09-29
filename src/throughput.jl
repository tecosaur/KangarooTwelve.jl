# This is a /little/ hacky(*), but this file is to be short-lived anyway.
# It's just for ad-hoc benchmarking.
#
# (*) Just a smiggin, the tinyiest bit 🤏. Shouldn't cause more
#  than a minor dent to your sanity if you care about small things
#  like not being a godawful hack.

function throughput(U::Type{<:Unsigned}, func, size, label=nameof(func))
    message = rand(U, round(Int, 1024*4^size ÷ sizeof(U)))
    start = time()
    x = func(message)
    stop = time()
    push!([], x) # prevent the call from being optimised away
    println(" $label throughput: ~$(round(Int, 1024*4^size / (stop - start) / 1024^2)) MiB/s")
end

throughput(func::Function, size, label=nameof(func)) =
    throughput(UInt8, func, size, label)

function throughput(::typeof(keccak_p1600), size)
    state = EMPTY_STATE
    rounds = round(Int, 1024*4^size / 200)
    start = time()
    for _ in 1:rounds
        state = keccak_p1600(state)
    end
    stop = time()
    push!([], state) # prevent the call from being optimised away
    println(" Keccak-p[1600,12] throughput: ~$(round(Int, 1024*4^size / (stop - start) / 1024^2)) MiB/s")
end

throughput(::typeof(turboshake), size) =
    throughput(UInt64, turboshake, size, "TurboSHAKE-128")

throughput(::typeof(k12_singlethreaded), size) =
    throughput(UInt64, k12_singlethreaded, size, "KangarooTwelve (singlethreaded)")

function throughput(::Val{:crc32c}, size)
    CRC32c_pkg = Base.PkgId(Base.UUID("8bf52ea8-c179-5cab-976a-9e18b702a9bc"), "CRC32c")
    CRC32c = Base.require(CRC32c_pkg)
    CRC32c.crc32c(UInt8[0x00]) # warmup :)
    throughput(CRC32c.crc32c, size)
end

function throughput(::Val{:sha1}, size)
    SHA_pkg = Base.PkgId(Base.UUID("ea8e919c-243c-51af-8825-aaa63cd721ce"), "SHA")
    SHA = Base.require(SHA_pkg)
    throughput(SHA.sha1, size)
end

function throughput(::Val{:sha256}, size)
    SHA_pkg = Base.PkgId(Base.UUID("ea8e919c-243c-51af-8825-aaa63cd721ce"), "SHA")
    SHA = Base.require(SHA_pkg)
    throughput(SHA.sha256, size)
end

function throughput(::Val{:sha3_256}, size)
    SHA_pkg = Base.PkgId(Base.UUID("ea8e919c-243c-51af-8825-aaa63cd721ce"), "SHA")
    SHA = Base.require(SHA_pkg)
    throughput(SHA.sha3_256, size)
end

function throughput(::Val{:md5}, size)
    MD5_pkg = Base.PkgId(Base.UUID("6ac74813-4b46-53a4-afec-0b5dc9d7885c"), "MD5")
    MD5 = Base.require(MD5_pkg)
    try
        MD5.md5(UInt8[0x00])
    catch e
        if e isa MethodError
            return Base.invokelatest(throughput, Val{:md5}(), size)
        else
            rethrow()
        end
    end
    throughput(MD5.md5, size)
end

throughput(alg) = throughput(alg, 9)

"""
    throughput(algorithm::Symbol, size=9)

Algorithms:
- keccak
- turboshake
- k12 (aka. k12_singlethreaded)
- crc32c
- sha1, sha256, sha3_256
- md5

|`size`| msg size |
|------|---------|
| 0    | 1 KiB   |
| 1    | 4 KiB   |
| 2    | 16 KiB  |
| 3    | 64 KiB  |
| 4    | 256 KiB |
| 5    | 1 MiB   |
| 5    | 4 MiB   |
| 7    | 16 MiB  |
| 8    | 64 MiB  |
| *9    | 256 MiB |
| 10   | 1 GiB   |
| 11 | 4 GiB   |
| ...  |         |
"""
function throughput(alg::Symbol, size=9) # 4 = 1KiB, 6 = 32KiB, 8 = 1MiB...
    func = if alg == :keccak
        keccak_p1600
    elseif alg == :turboshake
        turboshake
    elseif alg == :k12_singlethread
        k12_singlethreaded
    elseif alg == :k12
        k12
    elseif alg ∈ (:crc32, :crc32c)
        Val{:crc32c}()
    elseif alg ∈ (:sha1, :sha256, :sha3_256)
        Val(alg)
    elseif alg == :md5
        Val{:md5}()
    else
        error("I don't know that alg of the k12 algorithm.")
    end
    throughput(func, size)
end
