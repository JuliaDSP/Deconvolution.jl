
using NumberTheoreticTransforms
using LinearAlgebra
using SpecialMatrices

"""
    fermat(convolved, kernel)

Calculates deconvolution with Number Theoretic Transform.
"""
function fermat(convolved::AbstractArray{T, 1}, h::AbstractArray{T, 1}, g::T, q::T) where {T<:Integer}
    N = length(convolved)
    bias = div(q-1, 2) + 1 #negative numbers represeted by the upper half of [0, q) range
    H = fnt(h, g, q)
    H_inv = invmod.(H, q)
    hh = Circulant(h)
    D = det(hh) |> T #TODO: find determinant faster with FFT
    Dm = mod(D, q)
    ym = convolved
    YM = fnt(ym, g, q)
    X0 = mod.(H_inv .* YM, q)
    xm = ifnt(X0, g, q)
    x0 = mod.(Dm * xm, q)
    x0[x0 .>= bias] = x0[x0 .>= bias] .- q
    x = x0
    y_prim = hh * x0
    @assert rem.(ym * D - y_prim, q) == zeros(T, N)
    y = div.(ym * D - y_prim, q)
    ym = mod.(y, q)
    j = 1
    while y != zeros(T, length(convolved))
        YM = fnt(ym, g, q)
        Xj = H_inv .* YM
        xj = ifnt(Xj, g, q)
        xj[xj .>= bias] = xj[xj .>= bias] .- q
        x = x .+ (xj * q^j)
        y_prim = hh * xj
        @assert rem.(y - y_prim, q) == zeros(T, N)
        y = div.(y - y_prim, q)
        ym = mod.(y, q)
        j = j + 1
    end

    return x // D
end
