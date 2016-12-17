# Deconvolution.jl

| **Documentation**                       | [**Package Evaluator**][pkgeval-link] | **Build Status**                          | **Code Coverage**               |
|:---------------------------------------:|:-------------------------------------:|:-----------------------------------------:|:-------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![][pkg-0.5-img]][pkg-0.5-url]       | [![Build Status][travis-img]][travis-url] | [![][coveral-img]][coveral-url] |
| [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.6-img]][pkg-0.6-url]       | [![Build Status][appvey-img]][appvey-url] | [![][codecov-img]][codecov-url] |

Introduction
------------

This package provides a set of functions to
[deconvolve](https://en.wikipedia.org/wiki/Deconvolution) digital signals, like
images or time series.  This is written in [Julia](http://julialang.org/), a
modern high-level, high-performance dynamic programming language designed for
technical computing.

Installation
------------

`Deconvolution.jl` is available for Julia 0.5 and later versions, and can be
installed with
[Julia built-in package manager](http://docs.julialang.org/en/stable/manual/packages/).
In a Julia session run the command

```julia
julia> Pkg.update()
julia> Pkg.add("Deconvolution")
```

### Documentation

The complete manual of `Deconvolution.jl` is available at
http://deconvolutionjl.readthedocs.io.  It has more detailed explanation of the
methods used and the examples are complemented with pictures.  You can also
download the PDF version of the manual from
https://media.readthedocs.org/pdf/deconvolutionjl/latest/deconvolutionjl.pdf.

Usage
-----

Currently `Deconvolution.jl` provides only one methd, but others will come in
the future.

### `wiener`

```julia
wiener(input, signal, noise[, blurring])
```

The [Wiener deconvolution](https://en.wikipedia.org/wiki/Wiener_deconvolution)
attempts at reducing the noise in a digital signal by suppressing frequencies
with low
[signal-to-noise ratio](https://en.wikipedia.org/wiki/Signal-to-noise_ratio).
The signal is assumed to be degraded by additive noise and a shift-invariant
blurring function.

The `wiener` function can be used to apply the Wiener deconvolution method to a
digital signal.  The arguments are:

* `input`: the digital signal
* `signal`: the original signal (or a signal with a luckily similar power
  spectrum)
* `noise`: the noise of the signal (or a noise with a luckily similar power
  spectrum)
* `blurring` (optional argument): the blurring kernel

All arguments must be arrays, all with the same size, and all of them in the
time/space domain (they will be converted to the frequency domain internally
using `fft` function).  Argument `noise` can be also a real number, in which
case a constant noise with that value will be assumed (this is a good
approximation in the case of
[white noise](https://en.wikipedia.org/wiki/White_noise)).

Examples
--------

### Wiener deconvolution

Here is an example of use of `wiener` function to perform the Wiener
deconvolution of an image, degraded with a blurring function and an additive
noise.

``` julia
using Images, TestImages, Deconvolution, ImageView

# Open the test image
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
# Polish the image with Deconvolution deconvolution
polished = wiener(blurred_img, img, noise, blurring)

# Wiener deconvolution works also when you don't have the real image and noise,
# that is the most common and useful case.  This happens because the Wiener
# filter only cares about the power spectrum of the signal and the noise, so you
# don't need to have the exact signal and noise but something with a similar
# power spectrum.
img2 = float(data(testimage("livingroom"))) # Load another image
noise2 = rand(size(img)) # Create another additive noise
# Polish the image with Deconvolution deconvolution
polished2 = wiener(blurred_img, img2, noise2, blurring)

# Compare...
view(img) # ...the original image
view(blurred_img) # ...the blurred image
view(polished) # ...the polished image
view(polished2) # ...the second polished image
```

Development
-----------

The package is developed at https://github.com/JuliaDSP/Deconvolution.jl.  There
you can submit bug reports, propose new deconvolution methods with pull
requests, and make suggestions.  If you would like to take over maintainership
of the package in order to further improve it, please open an issue.

### History ###

The ChangeLog of the package is available in
[NEWS.md](https://github.com/JuliaDSP/Deconvolution.jl/blob/master/NEWS.md) file
in top directory.

License
-------

The `Deconvolution.jl` package is licensed under the MIT "Expat" License.  The
original author is Mosè Giordano.



[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://deconvolutionjl.readthedocs.io/en/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://deconvolutionjl.readthedocs.io/en/stable/

[pkgeval-link]: http://pkg.julialang.org/?pkg=Deconvolution

[pkg-0.5-img]: http://pkg.julialang.org/badges/Deconvolution_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/detail/Deconvolution.html
[pkg-0.6-img]: http://pkg.julialang.org/badges/Deconvolution_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/detail/Deconvolution.html

[travis-img]: https://travis-ci.org/JuliaDSP/Deconvolution.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaDSP/Deconvolution.jl

[appvey-img]: https://ci.appveyor.com/api/projects/status/8gfd4r6807w93umj/branch/master?svg=true
[appvey-url]: https://ci.appveyor.com/project/giordano/deconvolution-jl

[coveral-img]: https://coveralls.io/repos/github/JuliaDSP/Deconvolution.jl/badge.svg?branch=master
[coveral-url]: https://coveralls.io/github/JuliaDSP/Deconvolution.jl?branch=master

[codecov-img]: https://codecov.io/gh/JuliaDSP/Deconvolution.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaDSP/Deconvolution.jl
