### Deconvolution.jl ---  Deconvolution in Julia language
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
### Code:

module Deconvolution

using FFTW

include("wiener.jl")
include("lucy.jl")

end # module
