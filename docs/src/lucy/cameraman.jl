using Images, TestImages
using Deconvolution
using FFTW
using ZernikePolynomials #https://github.com/rdoelman/ZernikePolynomials.jl

img = channelview(testimage("cameraman"))

blurring = evaluateZernike(LinRange(-16,16,512), [12, 4, 0], [1.0, -1.0, 2.0], index=:OSA)
blurring = fftshift(blurring)
blurring = blurring ./ sum(blurring)

blurred_img = fft(img) .* fft(blurring) |> ifft |> real

@time result_200_iterations = lucy(blurred_img, blurring, iterations=200)
@time result_2000_iterations = lucy(blurred_img, blurring, iterations=2000)

(N, M) = size(blurred_img)

save("original.jpg", Images.clamp01nan.(img))
save("blurred.jpg", Images.clamp01nan.(blurred_img))
save("restored_200.jpg", Images.clamp01nan.(result_200_iterations))
save("restored_2000.jpg", Images.clamp01nan.(result_2000_iterations))
