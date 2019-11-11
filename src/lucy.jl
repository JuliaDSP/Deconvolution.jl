### lucy.jl ---  Lucy deconvolution
#
# Copyright (C) 2019 Jakub Wronowski.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: richadson, lucy, deconvolution, signal processing
#
# This file is a part of Deconvolution.jl.
#
# License is MIT "Expat".
#
### Commentary:
#
# This file provides the `lucy' function to perform Richardson-Lucy deconvolution.
#
### Code:

export lucy

function lucy(observed::AbstractArray, psf::AbstractArray; iterations::Integer=8)
    @assert size(observed) == size(psf)
    @assert iterations >= 0

    psf_ft = fft(psf)
    psf_ft_conj = conj.(psf_ft)

    function lucystep(e)
        return e .* real(ifft(fft(observed ./ ifft(fft(e) .* psf_ft)) .* psf_ft_conj)) 
    end

    estimated = observed
    for i in 1:iterations
        estimated = lucystep(estimated)
    end

    return estimated
end

"""
    lucy(observed, psf[, iterations])

Return the Richardson-Lucy deconvolution of `observed`, using the point spread
function `psf`.

All arguments must be arrays in the time/space domain and all of the same size,
they will be converted into the frequency domain internally using the `fft`
function.
"""
lucy
