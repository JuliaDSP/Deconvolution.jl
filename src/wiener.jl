### wiener.jl ---  Wiener deconvolution
#
# Copyright (C) 2016 Mosè Giordano.
#
# Maintainer: Mosè Giordano <mose AT gnu DOT org>
# Keywords: wiener, deconvolution, signal processing
#
# This file is a part of Deconvolution.jl.
#
# License is MIT "Expat".
#
### Commentary:
#
# This file provides the `wiener' function to perform Wiener deconvolution.  It
# has two modes: one for the case of presence of a known blurring, the other one
# without blurring.
#
### Code:

export wiener

### The implementation of the deconvolution algorithm

function _wiener_no_blur(Y::AbstractArray, # Fourier transform of the input
                         S::AbstractArray, # Power spectrum of the signal
                         N::AbstractArray) # Power spectrum of the noise
    # Without blurring, the filter is:
    #   |S|² / (|S|² + |N|²)
    return real(ifft(Y .* S ./ (S .+ N)))
end

function _wiener_with_blur(Y::AbstractArray, # Fourier transform of the input
                           S::AbstractArray, # Power spectrum of the signal
                           N::AbstractArray, # Power spectrum of the noise
                           H::AbstractArray) # Fourier transform of the blurring
    # With blurring, the filter is:
    #   H* / (|H|² + |N/S|²)
    return real(ifft(Y .* conj(H) ./ (abs2.(H) .+ N ./ S)))
end

### User interface

## Without blurring
function wiener(input::AbstractArray, signal::AbstractArray,
                noise::AbstractArray)
    @assert size(input) == size(signal) == size(noise)
    input_ft = fft(input)
    signal_power_spectrum = abs2.(fft(signal))
    noise_power_spectrum = abs2.(fft(noise))
    return _wiener_no_blur(input_ft, signal_power_spectrum,
                           noise_power_spectrum)
end

# Noise as a real (it's converted to an array)
function wiener(input::AbstractArray, signal::AbstractArray, noise::Real)
    @assert size(input) == size(signal)
    input_ft = fft(input)
    signal_power_spectrum = abs2.(fft(signal))
    noise_power_spectrum = noise .* ones(eltype(input), size(input))
    return _wiener_no_blur(input_ft, signal_power_spectrum,
                           noise_power_spectrum)
end

## With blurring
function wiener(input::AbstractArray, signal::AbstractArray,
                noise::AbstractArray, blurring::AbstractArray)
    @assert size(input) == size(signal) == size(noise) == size(blurring)
    input_ft = fft(input)
    signal_power_spectrum = abs2.(fft(signal))
    noise_power_spectrum = abs2.(fft(noise))
    blurring_ft = fft(blurring)
    return _wiener_with_blur(input_ft, signal_power_spectrum,
                             noise_power_spectrum, blurring_ft)
end

# Noise as a real (it's converted to an array)
function wiener(input::AbstractArray, signal::AbstractArray,
                noise::Real, blurring::AbstractArray)
    @assert size(input) == size(signal) == size(blurring)
    input_ft = fft(input)
    signal_power_spectrum = abs2.(fft(signal))
    noise_power_spectrum = noise .* ones(eltype(input), size(input))
    blurring_ft = fft(blurring)
    return _wiener_with_blur(input_ft, signal_power_spectrum,
                             noise_power_spectrum, blurring_ft)
end


"""
    wiener(input, signal, noise[, blurring])

Return the Wiener deconvolution of `input`, using the power spectrum of `signal`
and `noise`.  If the `input` was blurred with a known blurring kernel, pass it
as fourth argument, `blurring`.

All arguments must be arrays in the time/space domain and all of the same size,
they will be converted into the frequency domain internally using the `fft`
function.
"""
wiener
