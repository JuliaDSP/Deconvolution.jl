
using NumberTheoreticTransforms
using LinearAlgebra
using SpecialMatrices

"""
    fermat(convolved, kernel)

Calculates deconvolution with Number Theoretic Transform.
"""
function fermat(convolved::AbstractArray{T}, h::AbstractArray{T}, g::T, q::T) where {T<:Integer}
    H = fnt(h, g, q)
    H_inv = invmod.(H, q)
    hh = Circulant(h)
    @show D = det(hh) |> T #TODO: find determinant faster with FFT
    Dm = mod(D, q)
    ym0 = convolved
    YM = fnt(ym0, g, q)
    X0 = mod.(H_inv .* YM, q)
    xm = ifnt(X0, g, q)
    x0 = mod.(Dm * xm, q)
    x0[x0 .>= div(q-1, 2)+1] = x0[x0 .>= div(q-1, 2)+1] .- q
    x = x0
    y_prim = hh * x0
    ym0 * D - y_prim
    y = div.(ym0 * D - y_prim, q)
    @assert rem.(ym0 * D - y_prim, q) == zeros(T, length(y))
    ym = mod.(y, q)
    j = 1
    while y != zeros(T, length(convolved))
        YM = fnt(ym, g, q)
        Xj = H_inv .* YM
        xj = ifnt(Xj, g, q)
        xj[xj .>= div(q-1, 2)+1] = xj[xj .>= div(q-1, 2)+1] .- q
        x = x .+ (xj * q^j)
        y_prim = hh * xj
        y = div.(y - y_prim, q)
        ym = mod.(y, q)
        j = j + 1
    end

    return x // D
end
