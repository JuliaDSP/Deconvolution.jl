var documenterSearchIndex = {"docs":
[{"location":"#Deconvolution.jl","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"","category":"page"},{"location":"#Introduction","page":"Deconvolution.jl","title":"Introduction","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Deconvolution.jl provides a set of functions to deconvolve digital signals, like images or time series. This is written in Julia, a modern high-level, high-performance dynamic programming language designed for technical computing.","category":"page"},{"location":"#Installation","page":"Deconvolution.jl","title":"Installation","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The latest version of Deconvolution.jl is available for Julia 1.0 and later versions, and can be installed with Julia built-in package manager. In a Julia session, after entering the package manager mode with ], run the command","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"pkg> add Deconvolution","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Older versions are also available for Julia 0.4-0.7.","category":"page"},{"location":"#Usage","page":"Deconvolution.jl","title":"Usage","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Currently Deconvolution.jl provides only two methods, but others will hopefully come in the future.","category":"page"},{"location":"#wiener-function","page":"Deconvolution.jl","title":"wiener function","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The Wiener deconvolution attempts at reducing the noise in a digital signal by suppressing frequencies with low signal-to-noise ratio. The signal is assumed to be degraded by additive noise and a shift-invariant blurring function.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Theoretically, the Wiener deconvolution method requires the knowledge of the original signal, the blurring function, and the noise. However, these conditions are difficult to met (and, of course, if you know the original signal you do not need to perform a deconvolution in order to recover the signal itself), but a strength of the Wiener deconvolution is that it works in the frequency domain, so you only need to know with good precision the power spectra of the signal and the noise. In addition, most signals of the same class have fairly similar power spectra and the Wiener filter is insensitive to small variations in the original signal power spectrum. For these reasons, it is possible to estimate the original signal power spectrum using a representative of the class of signals being filtered.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"For a short review of the Wiener deconvolution method see http://www.dmf.unisalento.it/~giordano/allow_listing/wiener.pdf and references therein.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The wiener function can be used to apply the Wiener deconvolution method to a digital signal. The arguments are:","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"input: the digital signal\nsignal: the original signal (or a signal with a likely similar   power spectrum)\nnoise: the noise of the signal (or a noise with a likely similar   power spectrum)\nblurring (optional argument): the blurring kernel","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"All arguments must be arrays, all with the same size, and all of them in the time/space domain (they will be converted to the frequency domain internally using fft function). Argument noise can be also a real number, in which case a constant noise with that value will be assumed (this is a good approximation in the case of white noise).","category":"page"},{"location":"#lucy-function","page":"Deconvolution.jl","title":"lucy function","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The Richardson-Lucy deconvolution is an iterative method based on Bayesian inference for restoration of signal that is convolved with a point spread function.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The lucy function can be used to apply the Richardson-Lucy deconvolution method to a digital signal. The arguments are:","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"observed: the digital signal\npsf: the point spread function\niterations (optional argument): the number of iterations","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"First two arguments must be arrays, all with the same size, and all of them in the time/space domain (they will be converted to the frequency domain internally using fft function). Argument iterations is an integer number. The more iterations is specified the better result should be if the solution converges and it is going to converge if PSF is estimated well.","category":"page"},{"location":"#Examples","page":"Deconvolution.jl","title":"Examples","text":"","category":"section"},{"location":"#Wiener-deconvolution","page":"Deconvolution.jl","title":"Wiener deconvolution","text":"","category":"section"},{"location":"#Noisy-time-series","page":"Deconvolution.jl","title":"Noisy time series","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"This is an example of application of the Wiener deconvolution to a time series.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"We first construct the noisy signal:","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"using LombScargle, Deconvolution, Plots, Statistics\nt = range(0, stop=10, length=1000) # observation times\nx = sinpi.(t) .* cos.(5 .* t) - 1.5 .* cospi.(t) .* sin.(2 .* t) # the original signal\nn = rand(length(x)) # noise to be added\ny = x .+ 3 .* (n .- mean(n)) # observed noisy signal","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"In order to perform the Wiener deconvolution, we need a signal that has a power spectrum similar to that of the original signal. We can use the Lomb–Scargle periodogram to find out the dominant frequencies in the observed signal, as implemented in the the Julia package LombScargle.jl.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"# Lomb-Scargle periodogram\np = lombscargle(t, y, maximum_frequency=2, samples_per_peak=10)\nplot(freqpower(p)...)","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"After plotting the periodogram you notice that it has three peaks, one in each of the following intervals: 0 05, 05 1, 1 15. Use the LombScargle.model function to create the best-fitting Lomb–Scargle model at the three best frequencies, that can be found with the findmaxfreq function (see the manual at http://lombscarglejl.readthedocs.io/ for more details):","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"m1 = LombScargle.model(t, y, findmaxfreq(p, [0, 0.5])[1]) # first model\nm2 = LombScargle.model(t, y, findmaxfreq(p, [0.5, 1])[1]) # second model\nm3 = LombScargle.model(t, y, findmaxfreq(p, [1, 1.5])[1]) # third model","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Once you have these three frequencies, you can deconvolve y by feeding wiener with a simple signal that is the sum of these three models:","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"signal = m1 + m2 + m3 # signal for `wiener`\nnoise = rand(length(y)) # noise for `wiener`\npolished = wiener(y, signal, noise)\nplot(t, x, size=(900, 600), label=\"Original signal\", linewidth=2)\nplot!(t, y, label=\"Observed signal\")\nsavefig(\"time-series-observed.png\")\nplot(t, x, size=(900, 600), label=\"Original signal\", linewidth=2)\nplot!(t, polished, label=\"Recovered with Wiener\")\nplot!(t, signal, label=\"Lomb–Scargle model\")\nsavefig(\"time-series-recovered.png\")","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"(Image: image)","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"(Image: image)","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Note that the signal recovered with the Wiener deconvolution is generally a good improvement with respect to the best-fitting Lomb–Scargle model obtained using a few frequencies.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"With real-world data the Lomb–Scargle periodogram may not work as good as in this toy-example, but we showed a possible strategy to create a suitable signal to use with wiener function.","category":"page"},{"location":"#Blurred-noisy-time-series","page":"Deconvolution.jl","title":"Blurred noisy time series","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Additionally to noise, also a blurring kernel can applied to the image. First we define the blurring kernel as a Gaussian kernel. Using ifftshift we move the center to index position 1 which is required for the Wiener deconvolution algorithm.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"# Gaussian blurring kernel\nkernel = exp.( -10 .* (t .- 5) .^ 2)\nkernel ./= sum(kernel) # normalize kernel to sum of 1\nkernel = ifftshift(kernel) # move center to index pos 1","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The blurring can be applied simply in Fourier space. After the blurring we add the noise.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"y_blurred = real(ifft(fft(kernel) .* fft(x))) + noise","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The deblurred image can be obtained with the following call.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"deblurred = wiener(y_blurred, signal, noise, kernel)","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Additionally to the noise array we also pass kernel to wiener.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"plot(t, x, size=(900, 600), label=\"Original signal\", linewidth=2)\nplot!(t, deblurred, label=\"Deblurred with Wiener\")\nsavefig(\"time-series-deblurred.png\")","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"(Image: image)","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"(Image: image)","category":"page"},{"location":"#Blurred-image","page":"Deconvolution.jl","title":"Blurred image","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Here is an example of use of wiener function to perform the Wiener deconvolution of an image, degraded with a blurring function and an additive noise.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"using Images, TestImages, Deconvolution, ImageView\n\n# Open the test image\nimg = float(data(testimage(\"cameraman\")))'\n# Create the blurring kernel in frequency domain\nx = hcat(ntuple(x -> collect((1:512) - 257), 512)...)\nk = 0.001\nblurring_ft = exp.(-k .* (x .^ 2 .+ x' .^ 2) .^ (5//6))\n# Create additive noise\nnoise = rand(size(img))\n# Fourier transform of the blurred image, with additive noise\nblurred_img_ft = fftshift(blurring_ft) .* fft(img) + fft(noise)\n# Get the blurred image from its Fourier transform\nblurred_img = real(ifft(blurred_img_ft))\n# Get the blurring kernel in the space domain\nblurring = ifft(fftshift(blurring_ft))\n# Polish the image with Deconvolution deconvolution\npolished = wiener(blurred_img, img, noise, blurring)\n\n# Wiener deconvolution works also when you don't have the real image and noise,\n# that is the most common and useful case.  This happens because the Wiener\n# filter only cares about the power spectrum of the signal and the noise, so you\n# don't need to have the exact signal and noise but something with a similar\n# power spectrum.\nimg2 = float(data(testimage(\"livingroom\"))) # Load another image\nnoise2 = rand(size(img)) # Create another additive noise\n# Polish the image with Deconvolution deconvolution\npolished2 = wiener(blurred_img, img2, noise2, blurring)\n\n# Wiener also works using a real number instead of a noise array\npolished3 = wiener(blurred_img, img2, 1000, blurring)\npolished4 = wiener(blurred_img, img2, 10000, blurring)\n\n\n# Compare...\nview(img) # ...the original image\nview(blurred_img) # ...the blurred image\nview(polished) # ...the polished image\nview(polished2) # ...the second polished image\nview(polished3) # ...the third polished image\nview(polished4) # ...the fourth polished image","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Original image Blurred image\n(Image: ) (Image: )","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Image restored with exact power spectrum and noise Image restored with imperfect reference noise and spectrum\n(Image: ) (Image: )","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Image restored with imperfect spectrum and constant noise Image restored with imperfect spectrum and constant noise\n(Image: ) (Image: )","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Without knowing the noise array exactly the contrast drops significantly. Some postprocessing of the contrast can enhance the quality further.","category":"page"},{"location":"#Richardson-Lucy-deconvolution","page":"Deconvolution.jl","title":"Richardson-Lucy deconvolution","text":"","category":"section"},{"location":"#Blurred-image-2","page":"Deconvolution.jl","title":"Blurred image","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Here is an example of use of lucy function to perform the Richardson-Lucy deconvolution of an image convolved with point spread function that models lens aberration.","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"    using Images, TestImages, Deconvolution, FFTW, ZernikePolynomials, ImageView\n\n    img = channelview(testimage(\"cameraman\"))\n\n    # model of lens aberration\n    blurring = evaluateZernike(LinRange(-16,16,512), [12, 4, 0], [1.0, -1.0, 2.0], index=:OSA)\n    blurring = fftshift(blurring)\n    blurring = blurring ./ sum(blurring)\n\n    blurred_img = fft(img) .* fft(blurring) |> ifft |> real\n\n    @time restored_img_200 = lucy(blurred_img, blurring, iterations=200)\n    @time restored_img_2000 = lucy(blurred_img, blurring, iterations=2000)\n\n    imshow(img)\n    imshow(blurred_img)\n    imshow(restored_img_200)\n    imshow(restored_img_2000)","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Original image Blurred image\n(Image: ) (Image: )","category":"page"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"Result of 200 lucy iterations Result of 2000 lucy iterations\n(Image: ) (Image: )","category":"page"},{"location":"#Development","page":"Deconvolution.jl","title":"Development","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The package is developed at https://github.com/JuliaDSP/Deconvolution.jl. There you can submit bug reports, propose new deconvolution methods with pull requests, and make suggestions. If you would like to take over maintainership of the package in order to further improve it, please open an issue.","category":"page"},{"location":"#History","page":"Deconvolution.jl","title":"History","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The ChangeLog of the package is available in NEWS.md file in top directory.","category":"page"},{"location":"#License","page":"Deconvolution.jl","title":"License","text":"","category":"section"},{"location":"","page":"Deconvolution.jl","title":"Deconvolution.jl","text":"The Deconvolution.jl package is licensed under the MIT \"Expat\" License. The original author is Mosè Giordano.","category":"page"}]
}
