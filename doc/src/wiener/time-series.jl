### wiener-time-series.jl
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

using Random

Random.seed!(7)

using LombScargle, Deconvolution, Plots, Statistics, FFTW
t = range(0, stop=10, length=1000) # observation times
x = sinpi.(t) .* cos.(5t) - 1.5cospi.(t) .* sin.(2t) # the original signal

# Gaussian blurring kernel
kernel = exp.( - 10 .* (t .- 5).^2)  
kernel ./= sum(kernel) # normalize kernel to sum of 1
kernel = ifftshift(kernel) # move center to index pos 1

n = rand(length(x)) # noise to be added
noise = 3 .* (n .- mean(n))
y = x + noise # observed noisy signal
# blurred and noise signal
y_blurred = real(ifft(fft(kernel) .* fft(x))) + noise


# Lomb-Scargle periodogram
p = lombscargle(t, y, maximum_frequency=2, samples_per_peak=10)
plot(freqpower(p)...)

m1 = LombScargle.model(t, y, findmaxfreq(p, [0, 0.5])[1]) # first model
m2 = LombScargle.model(t, y, findmaxfreq(p, [0.5, 1])[1]) # second model
m3 = LombScargle.model(t, y, findmaxfreq(p, [1, 1.5])[1]) # third model

signal = m1 + m2 + m3
polished = wiener(y, signal, noise)
deblurred = wiener(y_blurred, signal,  noise, kernel)

# Plots 
plot(t, x, size=(900, 600), label="Original signal", linewidth=2)
plot!(t, y, label="Observed signal")
savefig("time-series-observed.png")

plot(t, x, size=(900, 600), label="Original signal", linewidth=2)
plot!(t, y_blurred, label="Blurred signal")
savefig("time-series-observed-blurred.png")

plot(t, x, size=(900, 600), label="Original signal", linewidth=2)
plot!(t, polished, label="Recovered with Wiener")
plot!(t, signal, label="Lomb–Scargle model")
savefig("time-series-recovered.png")

plot(t, x, size=(900, 600), label="Original signal", linewidth=2)
plot!(t, deblurred, label="Deblurred with Wiener")
savefig("time-series-deblurred.png")
