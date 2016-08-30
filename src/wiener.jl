### Wiener.jl ---  Wiener deconvolution
#
# Copyright (C) 2016 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: wiener, deconvolution, digital signals
#
# This file is a part of Wiener.jl.
#
# License is MIT "Expat".
#
### Code:

__precompile__()

module Wiener

export wiener

### The implementation of the algorithm

function _no_blur(Y::AbstractArray, # Fourier transform of the input
                  S::AbstractArray, # Power spectrum of the signal
                  N::AbstractArray) # Power spectrum of the noise
    # Without blurring, the filter is:
    #   |S|² / (|S|² + |N|²)
    return real(ifft(Y .* S ./ (S .+ N)))
end

function _with_blur(Y::AbstractArray, # Fourier transform of the input
                    S::AbstractArray, # Power spectrum of the signal
                    N::AbstractArray, # Power spectrum of the noise
                    H::AbstractArray) # Fourier transform of the blurring function
    # With blurring, the filter is:
    #   H* / (|H|² + |N/S|²)
    return real(ifft(Y .* conj(H) ./ (abs2(H) .+ N ./ S)))
end

### User interface

function wiener(input::AbstractArray, signal::AbstractArray,
                noise::AbstractArray)
    @assert size(input) == size(signal) == size(noise)
    input_ft = fft(input)
    signal_power_spectrum = abs2(fft(signal))
    noise_power_spectrum = abs2(fft(noise))
    return _no_blur(input_ft, signal_power_spectrum, noise_power_spectrum)
end

wiener(input::AbstractArray, signal::AbstractArray, noise::Real) =
    wiener(input, signal, ones(input)*noise)

function wiener(input::AbstractArray, signal::AbstractArray,
                noise::AbstractArray, blurring::AbstractArray)
    @assert size(input) == size(signal) == size(noise) == size(blurring)
    input_ft = fft(input)
    signal_power_spectrum = abs2(fft(signal))
    noise_power_spectrum = abs2(fft(noise))
    blurring_ft = fft(blurring)
    return _with_blur(input_ft, signal_power_spectrum,
                      noise_power_spectrum, blurring_ft)
end

wiener(input::AbstractArray, signal::AbstractArray,
       noise::Real, blurring::AbstractArray) =
           wiener(input, signal, ones(input)*noise, blurring)

"""
    wiener(input, signal, noise[, blurring])

Return the Wiener deconvolution of `input`, using the power spectrum of `signal`
and `noise`.  If the `input` was blurred with a known blurring function, pass it
as fourth argument, `blurring`.

All arguments must be in the time/space domain, they will be converted into the
frequency domain internally using the `fft` function.
"""
wiener

end # module
