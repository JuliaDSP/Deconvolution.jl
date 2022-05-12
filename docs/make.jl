using Documenter, Deconvolution

makedocs(modules = [Deconvolution],
         sitename = "Deconvolution.jl")

deploydocs(
    repo = "github.com/JuliaDSP/Deconvolution.jl.git",
    target = "build"
)
