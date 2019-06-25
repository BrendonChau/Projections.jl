module AffineTest

using Projections, Test, Random

@testset "affine" begin
Random.seed!(12345);
A = rand(5, 5)
b = ones(size(A, 1))
y = rand(5)
@test isapprox(
    project(Affine(A, b), y), 
    [0.329706,  -0.21199,   0.95313,   1.10365,  -0.36689],
    rtol = 1e-4)
end
end
