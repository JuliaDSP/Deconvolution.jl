using Images, TestImages
using Deconvolution
using FFTW
using ZernikePolynomials #https://github.com/rdoelman/ZernikePolynomials.jl

img = channelview(testimage("cameraman"))

blurring = evaluateZernike(LinRange(-16,16,512), [12, 4, 0], [1.0, -1.0, 2.0], index=:OSA)
blurring = blurring ./ maximum(blurring)
blurring = fftshift(blurring)
blurring = blurring ./ sum(blurring)

blurred_img = fft(img) .* fft(blurring) |> ifft |> real

@time result_200_iterations = lucy(blurred_img, blurring, iterations=200)
@time result_2000_iterations = lucy(blurred_img, blurring, iterations=2000)

(N, M) = size(blurred_img)

out = ((vcat(hcat(img, ones(M, 20), blurred_img),
             ones(20, 2M+20),
             hcat(result_200_iterations, ones(M, 20), result_2000_iterations))))

save("lucy-cameraman.jpg", Images.clamp01nan.(out))
