using Deconvolution
using Test, Random, FFTW

Random.seed!(42) # Fixed random seed

##### Wiener deconvolution

x = range(0, stop = 10, length = 15)
s = sinpi.(x) .- 1.5cos.(x)
n = rand(length(s))
@testset "No blurring" begin
    @test wiener(s + n, s, 0.5) ≈
        [-1.1111473161507488, -0.17041312679112874, -1.3469825775608049, 1.2447585407259756, 
          2.578860561637874, 1.09191362235022, 1.5339504088753226, -0.34582554585559355, 
          -1.316232389748769, -0.21271451029702804, -0.6157543266519075, -0.3338395657941805, 
          2.1195369797597614, 1.060270051101612, 1.4446941593151053]
end

@testset "Blurring" begin
    blurring_ft  = exp.(-0.001 * (range(-div(length(x), 2), length = length(x)) .^ 2) .^ (5 // 6))
    blurred_s_ft = fftshift(blurring_ft) .* fft(s) .+ fft(n)
    blurred_s    = real(ifft(blurred_s_ft))
    blurring     = ifft(fftshift(blurring_ft))
    @test wiener(blurred_s, s, 0.5, blurring) ≈
        [-1.1086114591808454, -0.17138676636801406, -1.3482585663111306, 1.2427663740234156, 
         2.5813314219359276, 1.0940885915582341, 1.532411498782748, -0.34834962964626176, 
         -1.311080379275804, -0.21869338550577228, -0.6084828281998079, -0.3366932366271657, 
         2.117735346098039, 1.0624233782615877, 1.4429704301669912]
end

##### Richardson-Lucy deconvolution

@testset "Richardson-Lucy deconvolution tests" begin
    k=-0.001
    blurring_ft  = exp.(-k * (range(-div(length(x), 2), length = length(x)) .^ 2) .^ (5 // 6))
    blurred_s_ft = fftshift(blurring_ft) .* fft(s)
    blurred_s    = real(ifft(blurred_s_ft))
    blurring     = ifft(fftshift(blurring_ft))

    estimated = lucy(blurred_s, blurring, iterations=1)

    @test sum(abs.(s .- estimated)) < sum(abs.(s .- blurred_s))

    k =-0.0015 # should improve also when blur params are not perfectly estimated 
    blurring_ft  = exp.(-k * (range(-div(length(x), 2), length = length(x)) .^ 2) .^ (5 // 6))
    blurring     = ifft(fftshift(blurring_ft))

    estimated = lucy(blurred_s, blurring, iterations=1)

    @test sum(abs.(s .- estimated)) < sum(abs.(s .- blurred_s))
end
