using KangarooTwelve
using Test

@testset "turboshake128" begin
    keccak_1600_init =
        [0xf1258f7940e1dde7, 0x84d5ccf933c0478a, 0xd598261ea65aa9ee, 0xbd1547306f80494d,
        0x8b284e056253d057, 0xff97a42d7f8e6fd4, 0x90fee5a0a44647c4, 0x8c5bda0cd6192e76,
        0xad30a6f71b19059c, 0x30935ab7d08ffc64, 0xeb5aa93f2317d635, 0xa9a6e6260d712103,
        0x81a57c16dbcf555f, 0x43b831cd0347c826, 0x01f22f1a11a5569f, 0x05e5635a21d9ae61,
        0x64befef28cc970f2, 0x613670957bc46611, 0xb87c5a554fd00ecb, 0x8c3ee88a1ccf32c8,
        0x940c7922ae3a2614, 0x1841f924a2c509e4, 0x16f53526e70465c2, 0x75f644e97f30a13b,
        0xeaf1ff7b5ceca249]

    keccak_1600_init12 =
        [0x32bdb29984dd0136, 0xd007ccf81390f00d, 0x255f7d637de64fab, 0x1e68c397bd8851b5,
        0x42650061d07bbd1d, 0x8904bc882d358cd1, 0x59861999977d5625, 0x3f58254b6d09b4b3,
        0xd16ca2220abb1d9c, 0x15f1da0b4c67b71b, 0xd65da28284ee7fb3, 0x286b5606107b116b,
        0xec032b5abf13af97, 0xb917b75d3d045350, 0x673d897d6f671820, 0x28d0476c209df0f0,
        0x0c6182b878ce4246, 0xb83f00cd9db7feff, 0x1a1419f9a37814c8, 0x29e53860a4af645f,
        0x6922db03f53e69e8, 0x36c30f12e8bb4606, 0x1a30b0200c06b2e2, 0xea5dc7c65491c36a,
        0xc1de6860d0d10358]

    keccak_1600_12rounds =
        [0x8b90f350b5018696, 0xbbd981435d1b7440, 0xc3907f0b1cf9fa9b, 0x650b8934855a881f,
        0x3acb12c9505d6125, 0x58aa8e6240715d20, 0x6dace6cf13225fcb, 0x17f2a9b889bad731,
        0x705552894792eb44, 0x39f3c319d2bfb4a0, 0x6ccbc0b5a1a26f58, 0xa59599376203c84f,
        0xb0a7f75e71b11faa, 0xf41131da8a5aceb4, 0xb780be5266849655, 0x67e882ed196882ca,
        0x669327075f0c94bb, 0x8a2c7727a5a3b612, 0x5e0d91393c42af39, 0xfa849e4b42de2f42,
        0x9c54fe0e1fd87d9b, 0x2bc0f29fd366b762, 0xe12aa82fc75e1704, 0x2d4bac302a4c445a,
        0x365f02c25ea643ee]

    keccak_1600_24rounds =
        [0x2d5c954df96ecb3c, 0x6a332cd07057b56d, 0x093d8d1270d76b6c, 0x8a20d9b25569d094,
        0x4f9c4f99e5e7f156, 0xf957b9a2da65fb38, 0x85773dae1275af0d, 0xfaf4f247c3d810f7,
        0x1f1b9ee6f79a8759, 0xe4fecc0fee98b425, 0x68ce61b6b9ce68a1, 0xdeea66c4ba8f974f,
        0x33c43d836eafb1f5, 0xe00654042719dbd9, 0x7cf8a9f009831265, 0xfd5449a6bf174743,
        0x97ddad33d8994b40, 0x48ead5fc5d0be774, 0xe3b8c8ee55b7b03c, 0x91a0226e649e42e9,
        0x900e3129e7badd7b, 0x202a9ec5faa3cce8, 0x5b3402464e1c3db6, 0x609f4e62a44c1059,
        0x20d06cd26a8fbf5c]

    @test KangarooTwelve.keccakP1600(copy(keccak_1600_init)) ==
        keccak_1600_12rounds
    @test KangarooTwelve.keccakP1600(copy(keccak_1600_init12)) ==
        keccak_1600_24rounds
    @test KangarooTwelve.keccakP1600(copy(keccak_1600_init), Val(24)) ==
        keccak_1600_24rounds
end
