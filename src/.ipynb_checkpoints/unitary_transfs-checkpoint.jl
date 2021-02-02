using RandomMatrices

function random_unitary_transf( matrix )
    dimension = length( matrix[:,1] )
    U = rand(Haar(2),dimension)
    U'*matrix*U
end

function infnts_random_unitary_transf( matrix, delta::Float64)
    dimension = length( matrix[:,1] )
    U1 = rand(Haar(2),dimension)
    exp_diag = Diagonal( exp.(1im*delta*(ones(dimension)-2*rand(dimension) ) )   )
    U = U1'*exp_diag*U1

    U'*matrix*U
end
