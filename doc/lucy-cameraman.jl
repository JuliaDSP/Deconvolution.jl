### lucy-cameraman.jl
#
# Copyright (C) 2019 Jakub Wronowski.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: deconvolution, signal processing
#
# This file is a part of Deconvolution.jl.
#
# License is MIT "Expat".
#
### Commentary:
#
# This file produces an example image used in Deconvolution.jl documentation.
#
### Code:

using Images, TestImages, Deconvolution, Random, FFTW

Random.seed!(42)

img = channelview(testimage("cameraman"))

(N, M) = size(img)
@assert N == M
@assert ispow2(N) # ensure size is power of 2

x = reshape(repeat(-N÷2:N÷2-1, M), (N,M))
k = 0.001

blurring_ft = @. exp(-k*(x ^ 2 + x ^ 2)^(5//6))
blurred_img_ft = fftshift(blurring_ft) .* fft(img)
blurred_img = ifft(blurred_img_ft) |> real
blurring = ifft(fftshift(blurring_ft)) |> real

exact_psf_result = lucy(blurred_img, blurring, iterations=100)

k=0.0007 # try a bit different blur
noise = rand(size(x)...)./20 # and add noise
blurred_img_ft = fftshift(blurring_ft) .* fft(img) .+ fft(noise)
blurred_img = ifft(blurred_img_ft) |> real
blurring_ft = @. exp(-k*(x ^ 2 + x ^ 2)^(5//6))
blurring = ifft(fftshift(blurring_ft)) |> real

different_psf_result = lucy(blurred_img, blurring, iterations=100)

out = ((vcat(hcat(img, ones(M, 20), blurred_img),
             ones(20, 2M+20),
             hcat(exact_psf_result, ones(M, 20), different_psf_result))))

save("lucy-cameraman.jpg", Images.clamp01nan.(out))



