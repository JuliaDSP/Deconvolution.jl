# Wiener

[![Build Status](https://travis-ci.org/giordano/Wiener.jl.svg?branch=master)](https://travis-ci.org/giordano/Wiener.jl)

[Wiener deconvolution](https://en.wikipedia.org/wiki/Wiener_deconvolution) in
[Julia language](http://julialang.org/).

``` julia
using Images, TestImages, Wiener, ImageView

# Open the image
img = float(data(testimage("cameraman")))'
# Create the blurring kernel in frequency domain
x = hcat(ntuple(x -> collect((1:512) - 257), 512)...)
k = 0.001
blurring_ft = exp(-k*(x .^ 2 + x' .^ 2).^(5//6))
# Create additive noise
noise = rand(size(img))
# Fourier transform of the blurred image, with additive noise
blurred_img_ft = fftshift(blurring_ft) .* fft(img) + fft(noise)
# Get the blurred image from its Fourier transform
blurred_img = real(ifft(blurred_img_ft))
# Get the blurring kernel in the space domain
blurring = ifft(fftshift(blurring_ft))
# Polish the image with Wiener deconvolution
polished = wiener(blurred_img, img, noise, blurring)

# Wiener deconvolution works also when you don't have the real image and noise,
# that is the most common and useful case.  This happens because the Wiener
# filter only cares about the power spectrum of the signal and the noise, so you
# don't need to have the exact signal and noise but something with a similar
# power spectrum.
img2 = float(data(testimage("livingroom"))) # Load another image
noise2 = rand(size(img)) # Create another additive noise
# Polish the image with Wiener deconvolution
polished2 = wiener(blurred_img, img2, noise2, blurring)

# Compare...
view(img) # ...the original image
view(blurred_img) # ...the blurred image
view(polished) # ...the polished image
view(polished2) # ...the second polished image
```
