@compile_workload begin
    k12(UInt8[])
    k12("keccak")
    k12(rand(UInt8, 8200))
    k12(rand(UInt8, 1024^2))
    k12(IOBuffer("keccak"))
    squeeze!([0x00], ByteSponge{21}())
end
