using Deconvolution
using Test, StableRNGs, FFTW

rng = StableRNG(42) # Fixed random seed

##### Wiener deconvolution

x = range(0, stop = 10, length = 15)
s = sinpi.(x) .- 1.5cos.(x)
n = rand(rng, length(s))
blurring_ft  = exp.(-0.001 * (range(-div(length(x), 2), length = length(x)) .^ 2) .^ (5 // 6))
blurred_s_ft = fftshift(blurring_ft) .* fft(s) .+ fft(n)
blurred_s    = real(ifft(blurred_s_ft))
blurring     = ifft(fftshift(blurring_ft))
@testset "No blurring with noise scalar parameter" begin
    @test wiener(s + n, s, 0.5) ≈
        [-1.0688616250172964, -0.32510758040709303, -0.4928689108095881, 1.644396339118468,
        1.866916504419297, 0.9245465020787589, 1.5461081346912358, 0.12677498594640266,
        -1.4157929359946475, -0.0247233772043046, -0.9345328207346199, 0.03628973013108201,
        2.414101393218902, 0.568743737331393, 1.5290604254094078]
end

@testset "Blurring with noise scalar parameter" begin
    @test wiener(blurred_s, s, 0.5, blurring) ≈
        [-1.0652058986210413, -0.3304166624246508, -0.4875855887501838, 1.6438301222753278,
        1.8625054116642619, 0.9303156760177144, 1.5413569689140112, 0.12849334093024348,
        -1.4141988631895932, -0.025180391629450346, -0.9343071603438413,
        0.03626923211984258, 2.4171155004554263, 0.5643293761857023, 1.5290637671248288]
end

@testset "No blurring with noise array" begin
    @test wiener(blurred_s, s, ifftshift(n)) ≈
        [-1.5060672134546482, -0.7416509948811038, -0.8868133826324668, 1.2108272657591024,
         1.414853150573584, 0.5580080918195482, 1.0351041993041612, -0.2750339384343115,
         -1.73737743125768, -0.5146485193142672, -1.3081993526690698, -0.37649294337342376,
         1.9491773714002598, 0.19435316342110073, 1.1019019675089878]
end

@testset "Blurring with noise array" begin
    @test wiener(blurred_s, s, ifftshift(n), blurring) ≈
        [-1.524937548601649, -0.7346444798363926, -0.8990054207462668, 1.2192518520141087,
        1.4207985817464395, 0.5541712011813319, 1.0415798812308625, -0.2742145640600878,
        -1.7491892004352156, -0.5071574425862486, -1.3165036261429737, -0.38315131571856315,
        1.9700562217921265, 0.18020610144859178, 1.1205661509755478]
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
