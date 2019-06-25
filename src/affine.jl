using LinearAlgebra

export Affine

struct Affine{T <: Real} <: ConvexSet
    VVt::Matrix{T}
    pinv_A_b::Vector{T}
    function Affine(
        A::Matrix{T}, 
        b::Vector{T}; 
        atol::Real = 0.0, 
        rtol::Real = (eps(real(float(one(T)))) * min(size(A)...)) * iszero(atol)
        ) where T <: Real
        size(A, 1) == length(b) || throw(ArgumentError(
            "Number of rows of A does not match length of b")
            )
        A = LinearAlgebra.svd!(A, full = false)
        tol = max(rtol * maximum(A.S), atol)
        Σ_pinv = zeros(length(A.S))
        index = A.S .> tol
        Σ_pinv[index] = 1. ./ A.S[index]
        Σ_pinv[findall(.!isfinite.(Σ_pinv))] .= 0.
        VVt = A.V * A.Vt
        display(Σ_pinv)
        pinv_A_b = A.V * Diagonal(Σ_pinv) * transpose(A.U) * b
        new{T}(VVt, pinv_A_b)
    end
end

function project!(s::Affine{T}, v::Vector{T}, y::Vector{T}) where T <: Real
    v .= y - s.VVt * y + pinv_A_b
end

function project(s::Affine{T}, y::Vector{T}) where T <: Real 
    project!(s, similar(y), y)
end