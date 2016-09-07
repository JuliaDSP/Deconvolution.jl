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

srand(7)

using LombScargle, Deconvolution, Plots
t = linspace(0, 10, 1000) # observation times
x = sinpi(t) .* cos(5t) - 1.5cospi(t) .* sin(2t) # the original signal
n = rand(length(x)) # noise to be added
y = x + 3(n - mean(n)) # observed noisy signal
# Lomb-Scargle periodogram
p = lombscargle(t, y, maximum_frequency=2, samples_per_peak=10)
plot(freqpower(p)...)

m1 = LombScargle.model(t, y, findmaxfreq(p, [0, 0.5])[1]) # first model
m2 = LombScargle.model(t, y, findmaxfreq(p, [0.5, 1])[1]) # second model
m3 = LombScargle.model(t, y, findmaxfreq(p, [1, 1.5])[1]) # third model

signal = m1 + m2 + m3
noise = rand(length(y)) # noise for `wiener`
polished = wiener(y, signal, noise)
plot(t, x, size=(900, 600), label="Original signal", linewidth=2)
plot!(t, y, label="Observed signal")
savefig("wiener-time-series-observed.png")
plot(t, x, size=(900, 600), label="Original signal", linewidth=2)
plot!(t, polished, label="Recovered with Wiener")
plot!(t, signal, label="Lomb–Scargle model")
savefig("wiener-time-series-recovered.png")
