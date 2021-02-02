using DensityMatrixModule

function partial_trace( density_matrix::DensityMatrix, qdit_index::DimType )

    new_num_qdits = density_matrix.num_qdits-1
    new_dims = zeros(Int8,new_num_qdits)

    skip=0
    for b=1:new_num_qdits
        if b==qdit_index
            skip = 1
            new_dims[b] = density_matrix.qdit_dimensions[b+skip]
            continue
        else
            new_dims[b] = density_matrix.qdit_dimensions[b+skip]
        end
    end

    new_prod_ds = prod_dims(new_dims)

    new_D = prod(new_dims[b] for b=1:new_num_qdits )
    new_rho = zeros(Complex{Float64}, new_D , new_D )


    for I=0:(new_D-1),J=0:(new_D-1)
        for i_qdit_index=0:(density_matrix.qdit_dimensions[qdit_index]-1)
            old_I = mod(I,density_matrix.prod_dimensions[qdit_index]) + i_qdit_index*density_matrix.prod_dimensions[qdit_index]
            old_J = mod(J,density_matrix.prod_dimensions[qdit_index]) + i_qdit_index*density_matrix.prod_dimensions[qdit_index]

            if qdit_index<=new_num_qdits
                old_I = old_I + sum( fld( mod(I,new_prod_ds[b+1] ) , new_prod_ds[b] )*density_matrix.prod_dimensions[b+1] 
                    for b=qdit_index:new_num_qdits )
                old_J = old_J + sum( fld( mod(J,new_prod_ds[b+1] ) , new_prod_ds[b] )*density_matrix.prod_dimensions[b+1] 
                    for b=qdit_index:new_num_qdits )
            end

            new_rho[I+1,J+1] = new_rho[I+1,J+1] + density_matrix.rho[old_I+1,old_J+1]
        end
    end

    construct_DensityMatrix(rho=new_rho,qdit_dimensions=new_dims,method=from_user)
end


function partial_trace(density_matrix ::DensityMatrix, qdit_indices ::Vector{DimType} )
    qdit_indices = sort(qdit_indices,rev=true)
    for qdit_index in qdit_indices
        density_matrix = partial_trace(density_matrix,qdit_index)
    end
    density_matrix
end
