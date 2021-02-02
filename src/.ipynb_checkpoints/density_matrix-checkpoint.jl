include("random_positive_matrix.jl")
include("unitary_transfs.jl")

#DimType=Int32
@enum DMConstructMethod from_user from_eigvals random_rank random_pure from_user_pure

struct DensityMatrix
    rho :: Matrix{Complex{Float64}}
    num_qdits :: Int64
    qdit_dimensions :: Vector{DimType}
    total_dimension :: DimType
    prod_dimensions :: Vector{DimType}
end

function construct_DensityMatrix(; qdit_dimensions::Vector{DimType}, 
            method::DMConstructMethod=from_user, 
            rho=nothing,
            rank=nothing,
            eigvals=nothing,
            U=nothing, 
            psi=nothing )
        
    #qdit_dimensions = qdit_dimensions
    num_qdits = length(qdit_dimensions)
    total_dimension = prod(qdit_dimensions[qdit] for qdit=1:num_qdits )
    prod_dimensions = [ 1 ; [ prod(qdit_dimensions[qdit_1] for qdit_1=1:qdit_2) for qdit_2=1:num_qdits ] ]

    if method==from_user
        rho = Hermitian( rho/tr(rho) )
    elseif method==from_eigvals
        if U==nothing
            rho = random_unitary_tranf( Diagonal(eigvals) )
        elseif U=="identity"
            rho = Diagonal(eigvals)
        else
            rho = U'*Diagonal(eigvals)*U
        end
    elseif method==random_pure
        rho = random_pure_rho( total_dimension )
    elseif method==from_user_pure
        rho = Hermitian((psi*psi')/(psi'*psi))
    elseif method==random_rank
        if rank==nothing
            rho = random_rho(total_dimension,total_dimension)
        else
            rho = random_rho(total_dimension,rank)
        end
    end
        
    return DensityMatrix( rho, num_qdits,qdit_dimensions,total_dimension, prod_dimensions )
end