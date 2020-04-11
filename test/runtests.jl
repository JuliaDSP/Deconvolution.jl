using Deconvolution
using Test, Random, FFTW

Random.seed!(42) # Fixed random seed

##### Wiener deconvolution

x = range(0, stop = 10, length = 15)
s = sinpi.(x) .- 1.5cos.(x)
n = rand(length(s))
blurring_ft  = exp.(-0.001 * (range(-div(length(x), 2), length = length(x)) .^ 2) .^ (5 // 6))
blurred_s_ft = fftshift(blurring_ft) .* fft(s) .+ fft(n)
blurred_s    = real(ifft(blurred_s_ft))
blurring     = ifft(fftshift(blurring_ft))
@testset "No blurring with noise scalar parameter" begin
    @test wiener(s + n, s, 0.5) ≈
        [-1.1111473161507488, -0.17041312679112874, -1.3469825775608049, 1.2447585407259756, 
          2.578860561637874, 1.09191362235022, 1.5339504088753226, -0.34582554585559355, 
          -1.316232389748769, -0.21271451029702804, -0.6157543266519075, -0.3338395657941805, 
          2.1195369797597614, 1.060270051101612, 1.4446941593151053]
end

@testset "Blurring with noise scalar parameter" begin
    @test wiener(blurred_s, s, 0.5, blurring) ≈
        [-1.1086114591808454, -0.17138676636801406, -1.3482585663111306, 1.2427663740234156, 
         2.5813314219359276, 1.0940885915582341, 1.532411498782748, -0.34834962964626176, 
         -1.311080379275804, -0.21869338550577228, -0.6084828281998079, -0.3366932366271657, 
         2.117735346098039, 1.0624233782615877, 1.4429704301669912]
end

@testset "No blurring with noise array" begin
    @test wiener(blurred_s, s, ifftshift(n)) ≈ [-1.4080797956919306, -0.49286583995619976, 
        -1.6205631520484833, 0.7973315264441291, 2.041487785980698, 0.7068033074341512, 
        1.303796832927927, -0.5691518774600562, -1.7799047917607576, -0.7841432287697161, 
        -0.9149839099122377, -0.4555342804458623, 1.761794603633166, 0.5975996848388563, 0.9542148100114622]
end

@testset "Blurring with noise array" begin
    @test wiener(blurred_s, s, ifftshift(n), blurring) ≈ [-1.4250670444336333, -0.48221282618655026, 
        -1.636981677852502, 0.8020690789792746, 2.053748123047759, 0.6974830808841243, 1.3177198257314031, 
        -0.5714565338646487, -1.7922287374907637, -0.7809316185676288, -0.9185450961060387, 
        -0.46118637814202684, 1.7785714975739182, 0.5887392863050276, 0.967947374977039]
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
