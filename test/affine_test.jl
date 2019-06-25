module AffineTest

using Projections, Test

@testset "affine" begin
A = reshape(collect(LinRange(0, 1, 25)), 5, 5)[:, 1:4]
b = collect(LinRange(1, 5, 5))
y = ones(size(A, 2))
@test isapprox(
    project(Affine(A, b), y), 
    [15.36,  9.12,  2.88,  -3.36],
    rtol = 1e-4)
end
end
