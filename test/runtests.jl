using KangarooTwelve
using Test

import KangarooTwelve: keccak_p1600, ingest, pad, squeeze, squeeze!,
    turboshake, Sponge, overwrite, ingest_length, k12_singlethreaded, k12

@testset "keccak" begin
    keccak_1600_init =
        (0xf1258f7940e1dde7, 0x84d5ccf933c0478a, 0xd598261ea65aa9ee, 0xbd1547306f80494d,
         0x8b284e056253d057, 0xff97a42d7f8e6fd4, 0x90fee5a0a44647c4, 0x8c5bda0cd6192e76,
         0xad30a6f71b19059c, 0x30935ab7d08ffc64, 0xeb5aa93f2317d635, 0xa9a6e6260d712103,
         0x81a57c16dbcf555f, 0x43b831cd0347c826, 0x01f22f1a11a5569f, 0x05e5635a21d9ae61,
         0x64befef28cc970f2, 0x613670957bc46611, 0xb87c5a554fd00ecb, 0x8c3ee88a1ccf32c8,
         0x940c7922ae3a2614, 0x1841f924a2c509e4, 0x16f53526e70465c2, 0x75f644e97f30a13b,
         0xeaf1ff7b5ceca249)

    keccak_1600_init12 =
        (0x32bdb29984dd0136, 0xd007ccf81390f00d, 0x255f7d637de64fab, 0x1e68c397bd8851b5,
         0x42650061d07bbd1d, 0x8904bc882d358cd1, 0x59861999977d5625, 0x3f58254b6d09b4b3,
         0xd16ca2220abb1d9c, 0x15f1da0b4c67b71b, 0xd65da28284ee7fb3, 0x286b5606107b116b,
         0xec032b5abf13af97, 0xb917b75d3d045350, 0x673d897d6f671820, 0x28d0476c209df0f0,
         0x0c6182b878ce4246, 0xb83f00cd9db7feff, 0x1a1419f9a37814c8, 0x29e53860a4af645f,
         0x6922db03f53e69e8, 0x36c30f12e8bb4606, 0x1a30b0200c06b2e2, 0xea5dc7c65491c36a,
         0xc1de6860d0d10358)

    keccak_1600_12rounds =
        (0x8b90f350b5018696, 0xbbd981435d1b7440, 0xc3907f0b1cf9fa9b, 0x650b8934855a881f,
         0x3acb12c9505d6125, 0x58aa8e6240715d20, 0x6dace6cf13225fcb, 0x17f2a9b889bad731,
         0x705552894792eb44, 0x39f3c319d2bfb4a0, 0x6ccbc0b5a1a26f58, 0xa59599376203c84f,
         0xb0a7f75e71b11faa, 0xf41131da8a5aceb4, 0xb780be5266849655, 0x67e882ed196882ca,
         0x669327075f0c94bb, 0x8a2c7727a5a3b612, 0x5e0d91393c42af39, 0xfa849e4b42de2f42,
         0x9c54fe0e1fd87d9b, 0x2bc0f29fd366b762, 0xe12aa82fc75e1704, 0x2d4bac302a4c445a,
         0x365f02c25ea643ee)

    keccak_1600_24rounds =
        (0x2d5c954df96ecb3c, 0x6a332cd07057b56d, 0x093d8d1270d76b6c, 0x8a20d9b25569d094,
         0x4f9c4f99e5e7f156, 0xf957b9a2da65fb38, 0x85773dae1275af0d, 0xfaf4f247c3d810f7,
         0x1f1b9ee6f79a8759, 0xe4fecc0fee98b425, 0x68ce61b6b9ce68a1, 0xdeea66c4ba8f974f,
         0x33c43d836eafb1f5, 0xe00654042719dbd9, 0x7cf8a9f009831265, 0xfd5449a6bf174743,
         0x97ddad33d8994b40, 0x48ead5fc5d0be774, 0xe3b8c8ee55b7b03c, 0x91a0226e649e42e9,
         0x900e3129e7badd7b, 0x202a9ec5faa3cce8, 0x5b3402464e1c3db6, 0x609f4e62a44c1059,
         0x20d06cd26a8fbf5c)

    @test keccak_p1600(keccak_1600_init) == keccak_1600_12rounds
    @test keccak_p1600(keccak_1600_init12) == keccak_1600_24rounds
    @test keccak_p1600(keccak_1600_init, Val(24)) == keccak_1600_24rounds
end

@testset "Sponge ingestion" begin
    sponge = ingest(Sponge(), 0x11)
    rate = (((::Sponge{rate}) where {rate}) -> rate)(sponge)
    @test sponge.state[1] == 0x0000000000000011
    sponge = ingest(sponge, 0x2222)
    @test sponge.state[1] == 0x0000000000222211
    sponge = ingest(sponge, 0x33333333)
    @test sponge.state[1] == 0x0033333333222211
    sponge = ingest(sponge, 0x4444)
    @test sponge.state[1] == 0x4433333333222211
    @test sponge.state[2] == 0x0000000000000044
    sponge = ingest(sponge, 0x5555555555555555)
    @test sponge.state[2] == 0x5555555555555544
    @test sponge.state[3] == 0x0000000000000055
    sponge = ingest(sponge, 0x66, 0x6666, 0x66666666)
    @test sponge.state[3] == 0x6666666666666655
    @test sponge.state[4] == 0x0000000000000000
    for _ in 4:rate-1
        sponge = ingest(sponge, 0x1234567812345678)
    end
    @test sum(sponge.state) == 0x3568ace824579ba2
    @test ingest(sponge, 0x1111111111111111).state[1] == 0x070513f3bdbfaa6f
    @test ingest(sponge, 0x22222222, 0x1111111111111111).state[1] == 0x0000000011111111
end

@testset "Sponge squeezing" begin
    sponge = Sponge{21}(reinterpret(NTuple{25, UInt64}, Tuple((UInt8(i) for i in 1:200))), 1)
    @test squeeze(UInt8, sponge) == 0x01
    @test squeeze(UInt16, sponge) == 0x0201
    @test squeeze(UInt32, sponge) == 0x04030201
    @test squeeze(UInt64, sponge) == 0x0807060504030201
    @test squeeze(UInt128, sponge) == 0x100f0e0d0c0b0a090807060504030201
    # Representational equivalence
    @test reinterpret(UInt128, squeeze(NTuple{2, UInt64}, sponge)) == 0x100f0e0d0c0b0a090807060504030201
    @test reinterpret(UInt128, squeeze(NTuple{4, UInt32}, sponge)) == 0x100f0e0d0c0b0a090807060504030201
    @test reinterpret(UInt128, squeeze(NTuple{8, UInt16}, sponge)) == 0x100f0e0d0c0b0a090807060504030201
    @test reinterpret(UInt128, squeeze(NTuple{16, UInt8}, sponge)) == 0x100f0e0d0c0b0a090807060504030201
    # Wrap around
    @test squeeze(NTuple{22, UInt64}, sponge)[end] == first(keccak_p1600(sponge.state))
    @test squeeze!(Vector{UInt64}(undef, 1000), sponge)[22] == first(keccak_p1600(sponge.state))
end

@testset "Length encoding" begin
    @test ingest_length(Sponge(), 0).state[1] == 0x0000000000000000
    @test ingest_length(Sponge(), 12).state[1] == 0x000000000000010c
    @test ingest_length(Sponge(), 65538).state[1] == 0x0000000003020001
end

bitpattern(num::Int) =
    Iterators.take(Iterators.cycle(0x00:0xfa), num) |> collect

@testset "Turboshake" begin
    # Generated from the reference implementation: TurboSHAKE.py
    # 256-bit capacity
    @test turboshake(UInt8,   UInt8[],          0x01, Val(256)) == 0x86
    @test turboshake(UInt16,  UInt8[],          0x01, Val(256)) == 0x8c86
    @test turboshake(UInt32,  UInt8[],          0x01, Val(256)) == 0x53bd8c86
    @test turboshake(UInt64,  UInt8[],          0x01, Val(256)) == 0x5a2078b053bd8c86
    @test turboshake(UInt128, UInt8[],          0x01, Val(256)) == 0x037d1f945d8185bb5a2078b053bd8c86
    @test turboshake(UInt128, UInt8[],          0x07, Val(256)) == 0x0f43edfc8c0443a2668c3b0bd33a225a
    @test turboshake(UInt128, UInt8[],          0x0c, Val(256)) == 0x3b3a8b4dae919b98fad5d126e862642c
    @test turboshake(UInt128, UInt8[],          0x17, Val(256)) == 0xfcb37fc2b8b697475bac3dcdebd5bf69
    @test turboshake(UInt128, UInt8[],          0x80, Val(256)) == 0x8237b3b479722729e95a8b6d8463efde
    @test turboshake(UInt128, bitpattern(17^1), 0x01, Val(256)) == 0xbaa8f812d0975b34ed14710a335f0f6f
    @test turboshake(UInt128, bitpattern(17^2), 0x01, Val(256)) == 0x819297eb5b6ee1b0adb55373a3ca3262
    @test turboshake(UInt128, bitpattern(17^3), 0x01, Val(256)) == 0x68a03a5687847180aae2860787058166
    @test turboshake(UInt128, bitpattern(17^4), 0x01, Val(256)) == 0xd6acc24a78d14501c296c50edde75d79
    @test turboshake(UInt128, bitpattern(17^5), 0x01, Val(256)) == 0x9107718af0504ff7f7bcbc6252e08541
    @test turboshake(UInt128, bitpattern(17^6), 0x01, Val(256)) == 0xef49a527ffd73f43aa8ff5f76c0ce4db
    # 512-bit capacity
    @test turboshake(UInt8,   UInt8[],          0x01, Val(512)) == 0xe3
    @test turboshake(UInt16,  UInt8[],          0x01, Val(512)) == 0xdde3
    @test turboshake(UInt32,  UInt8[],          0x01, Val(512)) == 0xf02ddde3
    @test turboshake(UInt64,  UInt8[],          0x01, Val(512)) == 0x6dde3b94f02ddde3
    @test turboshake(UInt128, UInt8[],          0x01, Val(512)) == 0x5cf35960c39ee3826dde3b94f02ddde3
    @test turboshake(UInt128, UInt8[],          0x07, Val(512)) == 0x49d0d015955ccf8c53f1f8ec065b554a
    @test turboshake(UInt128, UInt8[],          0x0c, Val(512)) == 0x16f64c66851915a60695f15745a8783c
    @test turboshake(UInt128, UInt8[],          0x17, Val(512)) == 0xba1d480c64e91d71f11ff3320242f6dc
    @test turboshake(UInt128, UInt8[],          0x80, Val(512)) == 0x3f8f58909ad1754390bbc184dfa4032a
    @test turboshake(UInt128, bitpattern(17^1), 0x01, Val(512)) == 0xbd5b67f2a842a20753b75587187da41d
    @test turboshake(UInt128, bitpattern(17^2), 0x01, Val(512)) == 0x909b27e2294e769db016f97087938ca4
    @test turboshake(UInt128, bitpattern(17^3), 0x01, Val(512)) == 0xc2c23fd3c73a5cc7a7ba463a8d66e875
    @test turboshake(UInt128, bitpattern(17^4), 0x01, Val(512)) == 0x01e39bd978c2113fa37b0ce45396a4ff
    @test turboshake(UInt128, bitpattern(17^5), 0x01, Val(512)) == 0x112daf7f8f5e9dfa401867b8beb3d22a
    @test turboshake(UInt128, bitpattern(17^6), 0x01, Val(512)) == 0x9843bab8cbb9f6ee85954968a7b32c61
    # Odd capacities
    # (I suspect these are broken because the capacities aren't a multiple of 64)
    @test_broken turboshake(UInt128, bitpattern(11^4), 0x2a, Val(8))   == 0x078cfb89b5b84ba84d30b26a92e48740
    @test_broken turboshake(UInt128, bitpattern(11^4), 0x22, Val(56))  == 0xfd731d4a42658ff64f327ae0a6f76fef
    @test_broken turboshake(UInt128, bitpattern(11^4), 0x7e, Val(520)) == 0x86a385f5f61243bc0b8f6a9f36e874d4
    # Representational equivalence
    let data = [UInt16(n) for n in 1:10000]
        @test turboshake(UInt128, reinterpret(UInt8, data))  == 0xede558468dc82194f78e10a5658c8309
        @test turboshake(UInt128, reinterpret(UInt16, data)) == 0xede558468dc82194f78e10a5658c8309
        @test turboshake(UInt128, reinterpret(UInt32, data)) == 0xede558468dc82194f78e10a5658c8309
        @test turboshake(UInt128, reinterpret(UInt64, data)) == 0xede558468dc82194f78e10a5658c8309
    end
end

@testset "Turboshake SIMD" begin
    # As before, just now x4
    x4(x) = (x, x, x, x)
    bitpattern64(n) = reinterpret(UInt64, bitpattern(8 * n))
    @test turboshake(UInt128, x4(UInt64[]),           0x01, Val(256)) == x4(turboshake(UInt128, UInt64[],           0x01, Val(256)))
    @test turboshake(UInt128, x4(bitpattern64(15)),   0x01, Val(256)) == x4(turboshake(UInt128, bitpattern64(15),   0x01, Val(256)))
    @test turboshake(UInt128, x4(bitpattern64(15^2)), 0x01, Val(256)) == x4(turboshake(UInt128, bitpattern64(15^2), 0x01, Val(256)))
    @test turboshake(UInt128, x4(bitpattern64(15^3)), 0x01, Val(256)) == x4(turboshake(UInt128, bitpattern64(15^3), 0x01, Val(256)))
    @test turboshake(UInt128, x4(bitpattern64(15^4)), 0x01, Val(256)) == x4(turboshake(UInt128, bitpattern64(15^4), 0x01, Val(256)))
    @test turboshake(UInt128, x4(bitpattern64(15^5)), 0x01, Val(256)) == x4(turboshake(UInt128, bitpattern64(15^5), 0x01, Val(256)))
    @test turboshake(UInt128, x4(bitpattern64(15^6)), 0x01, Val(256)) == x4(turboshake(UInt128, bitpattern64(15^6), 0x01, Val(256)))
end
