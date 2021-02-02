using RandomMatrices
using LinearAlgebra

DimType=Int

function random_pure_rho(dimension::DimType)
    psi = randn(dimension) + randn(dimension)*1im
    rho = (psi*psi')/(psi'*psi)
    Hermitian(rho)
end

function random_rho(dimension::DimType,rank::DimType)
    U = rand(Haar(2),dimension)
    Apos = U'*Diagonal( [rand(rank);zeros(dimension-rank)] )*U

    Hermitian(Apos/tr(Apos))
end
