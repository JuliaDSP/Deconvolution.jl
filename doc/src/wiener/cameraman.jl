### wiener-cameraman.jl
#
# Copyright (C) 2016 Mosè Giordano.
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

# Open the test image
img = channelview(testimage("cameraman"))
# Create the blurring kernel in frequency domain
x = hcat(ntuple(x -> collect((1:512) .- 257), 512)...)
k = 0.001
blurring_ft = @. exp(-k*(x ^ 2 + x ^ 2)^(5//6))
# Create additive noise
noise = rand(Float64, size(img))
# Fourier transform of the blurred image, with additive noise
blurred_img_ft = fftshift(blurring_ft) .* fft(img) .+ fft(noise)
# Get the blurred image from its Fourier transform
blurred_img = real(ifft(blurred_img_ft))
# Get the blurring kernel in the space domain
blurring = ifft(fftshift(blurring_ft))
# Polish the image with Deconvolution deconvolution
polished = wiener(blurred_img, img, noise, blurring)

# Wiener deconvolution works also when you don't have the real image and noise,
# that is the most common and useful case.  This happens because the Wiener
# filter only cares about the power spectrum of the signal and the noise, so you
# don't need to have the exact signal and noise but something with a similar
# power spectrum.
img2 = channelview(testimage("livingroom")) # Load another image
noise2 = rand(Float64, size(img)) # Create another additive noise
# Polish the image with Deconvolution deconvolution
polished2 = wiener(blurred_img, img2, noise2, blurring)

save("original.jpg", Images.clamp01nan.(img))
save("blurred.jpg", Images.clamp01nan.(blurred_img))
save("polished.jpg", Images.clamp01nan.(polished))
save("polished2.jpg", Images.clamp01nan.(polished2))
