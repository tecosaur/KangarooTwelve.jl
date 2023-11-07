#!/usr/bin/env -S julia
using SHA
using CRC32c
using MD5
using Blake3Hash
using KangarooTwelve: k12_singlethreaded, k12_singlethreaded_simd, k12_multithreaded, BLOCK_SIZE
using .GC

using CairoMakie
using ColorSchemes

function humansize(bytes::Integer; digits::Int=1)
    units = ("B", "KiB", "MiB", "GiB", "TiB", "PiB")
    magnitude = floor(Int, log(1024, max(1, bytes)))
    if 1024 <= bytes < 10.0^(digits-1) * 1024^magnitude
        magdigits = floor(Int, log10(bytes / 1024^magnitude)) + 1
        round(bytes / 1024^magnitude; digits = digits - magdigits)
    else
        round(Int, bytes / 1024^magnitude)
    end, units[1+magnitude]
end

bitpattern(num::Int) = Iterators.take(Iterators.cycle(0x00:0xfa), num) |> collect

alg_hashers = Dict(
    :sha256 => sha256,
    :sha3_256 => sha3_256,
    :md5 => md5,
    :crc32c => crc32c,
    :k12_singlethreaded => Base.Fix2(k12_singlethreaded, UInt[]),
    :k12_singlethreaded_simd => Base.Fix2(k12_singlethreaded_simd, UInt[]),
    :k12_multithreaded => Base.Fix2(k12_multithreaded, UInt[]),
    :blake3 => function (data)
        hasher = Blake3Ctx()
        Blake3Hash.update!(hasher, data)
        digest(hasher)
    end
)

function bench(alg::Symbol, size::Int; repeats::Int=clamp(round(Int, log2(2^30/size))^3, 5, 10000))
    alg == :k12_multithreaded && size <= BLOCK_SIZE && return NaN
    alg == :k12_singlethreaded_simd && size <= 4 * BLOCK_SIZE && return NaN
    hashfn = alg_hashers[alg]
    hashfn([0x01]) # Potentially trigger JIT
    data = bitpattern(size)
    batchsize = max(1, round(Int, sqrt(repeats ÷ 20)))
    times = UInt64[]
    GC.enable(false)
    for _ in 1:batchsize:repeats
        start = time_ns()
        for _ in 1:batchsize
            x = hashfn(data)
            Base.donotdelete(x)
        end
        stop = time_ns()
        push!(times, stop - start)
    end
    empty!(data)
    GC.enable(true)
    GC.gc()
    minimum(times) / (batchsize * 10^9)
end

all_algs = [:crc32c, :sha256, :sha3_256, :md5, :blake3, :k12_singlethreaded, :k12_multithreaded] # add SIMD once it's better
data_sizes = 2 .^ (10:30)

alg_labels = Dict(:k12_singlethreaded => "k12 (singlethreaded)",
                  :k12_singlethreaded_simd => "k12 (singlethreaded, SIMDx4)",
                  :k12_multithreaded => "k12 (multithreaded)",
                  :sha256 => "SHA2-256", :sha3_256 => "SHA3-256")

alg_colors =  Dict(
    :crc32c             => ColorSchemes.tol_bright.colors[end],
    :sha256             => ColorSchemes.tol_rainbow.colors[21],
    :sha3_256           => ColorSchemes.tol_rainbow.colors[24],
    :md5                => ColorSchemes.tol_rainbow.colors[7],
    :blake3             => ColorSchemes.tol_rainbow.colors[12],
    :k12_singlethreaded => ColorSchemes.tol_rainbow.colors[16],
    :k12_singlethreaded_simd => ColorSchemes.tol_rainbow.colors[17],
    :k12_multithreaded  => ColorSchemes.tol_rainbow.colors[15],
)

function scalingdata(algs::Vector{Symbol} = all_algs, dsizes::Vector{Int} = data_sizes)
    timings = Float64[]
    sizes = Int[]
    algorithms = Symbol[]
    for alg in algs
        @info "Benchmarking $alg"
        for size in dsizes
            time = bench(alg, size)
            push!(timings, time)
            push!(sizes, size)
            push!(algorithms, alg)
        end
    end
    (; timings, sizes, algorithms)
end

function scalingplot(data::NamedTuple)
    fig = Figure(resolution = (1200, 800))
    ax = Axis(fig[1,1], xscale = log2, yscale = log10,
              xlabel = "Input size", ylabel = "Throughput (GiB/s)",
              xtickformat = xs -> join.(humansize.(Int.(xs))),
              yticks = [0.1, 0.3, 1.0, 3.0, 10.0, 30.0],
              ytickformat = ys -> map(y -> if y < 1 string(y) else string(round(Int, y)) end, ys))
    algdata = map(all_algs ∩ unique(data.algorithms)) do alg
        mask = data.algorithms .== alg
        timings = data.timings[mask]
        sizes = data.sizes[mask]
        alg => (; timings, sizes)
    end
    sort!(algdata, by = ((_, data),) -> - last(data.sizes) / last(data.timings))
    for (alg, data) in algdata
        scatterlines!(ax, data.sizes, data.sizes ./ data.timings ./ 1024^3,
                      label = String(get(alg_labels, alg, alg)),
                      color = get(alg_colors, alg, ColorSchemes.tol_bright.colors[end]))
    end
    fig[1, 2] = Legend(fig, ax, "Algorithm", framevisible = false)
    fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    data = scalingdata()
    @info "Plotting"
    Makie.set_theme!(
        textcolor = :gray45,
        linecolor = :gray60,
        backgroundcolor=:transparent,
        Axis=(backgroundcolor=:transparent,
              xgridcolor = (:grey60, 0.3),
              ygridcolor = (:grey60, 0.3),
              leftspinevisible = false,
              rightspinevisible = false,
              bottomspinevisible = false,
              topspinevisible = false,
              xticksvisible = false,
              yticksvisible = false,),
        Legend=(bgcolor=:transparent,
                framevisible=false))
    save("scaling-benchmark.svg", scalingplot(data))
end
