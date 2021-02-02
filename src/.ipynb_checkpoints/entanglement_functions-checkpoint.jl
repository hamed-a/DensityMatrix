using LinearAlgebra
using DensityMatrixModule

const EPSIL = 1e-10

function entropy(rho)
    evals = real( eigvals(rho) )
    entropy = -sum( ( ( evals[i]<EPSIL ) ? 0 :  evals[i]*log(2, evals[i] )  ) for i in 1:length(evals) )
    entropy
end


function log_negativity(density_matrix::DensityMatrix,qdit_indices )
    transposed_density_matrix = partial_transpose( density_matrix,qdit_indices )
    evals = real( eigvals(transposed_density_matrix.rho) )
    log(2,sum(abs.( evals )  )  )
end


function negativity(density_matrix::DensityMatrix, qdit_indices )
    transposed_density_matrix = partial_transpose( density_matrix,qdit_indices )
    evals = real( eigvals(transposed_density_matrix.rho) )
    (sum(abs.( evals )  )-1)/2
end
