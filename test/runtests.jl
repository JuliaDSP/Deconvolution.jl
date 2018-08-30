using Deconvolution
using Test, Random, FFTW

Random.seed!(42) # Fixed random seed

##### Wiener deconvolution

x = range(0, stop = 10, length = 15)
s = sinpi.(x) .- 1.5cos.(x)
n = rand(length(s))

@testset "No blurring" begin
    @test wiener(s + n, s, 0.5) ≈
        [-1.5300909616077885,-0.4607553796475224,-1.73313390900302,0.8555916473930066,
         2.2692599396286277,0.7989152977213652,1.1431959877178226,-0.8118586562136341,
         -1.6519310564150556,-0.4882369786466372,-1.0387206070662094,
         -0.8324097927723406,1.8208048770862535,0.7427114498907227,1.0587917379150844]
end

@testset "Blurring" begin
    blurring_ft  = exp.(-0.001 * (range(-div(length(x), 2), length = length(x)) .^ 2) .^ (5 // 6))
    blurred_s_ft = fftshift(blurring_ft) .* fft(s) .+ fft(n)
    blurred_s    = real(ifft(blurred_s_ft))
    blurring     = ifft(fftshift(blurring_ft))
    @test wiener(blurred_s, s, 0.5, blurring) ≈
        [-1.5288865581538547,-0.4592713125163481,-1.7360819833552192,0.8525544925861452,
         2.272966823544607,0.8021732396525801,1.14091217864045,-0.8166101187683248,
         -1.6455701546601862,-0.4926325667795202,-1.0323127293402468,
         -0.8380228416032021,1.8214493186523952,0.7447075584686396,1.0566048110361044]
end
