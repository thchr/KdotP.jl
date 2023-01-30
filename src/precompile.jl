using SnoopPrecompile

@precompile_setup begin
    lgir = realify(lgirreps(230)["P"])[1]
    io = IOBuffer()
    @precompile_all_calls begin

        H = kdotp(lgir; timereversal=true, degree=3)
        show(io, MIME"text/plain"(), H)

        H[2]([.1,.1,.1],1)
    end
end