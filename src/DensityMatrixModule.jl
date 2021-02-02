module DensityMatrixModule
include("density_matrix.jl")
export DensityMatrix, DimType, DMConstructMethod, construct_DensityMatrix
export from_user, from_eigvals, random_rank, random_pure, from_user_pure

end
