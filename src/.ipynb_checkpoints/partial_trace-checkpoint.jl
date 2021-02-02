using DensityMatrixModule

function partial_trace( density_matrix::DensityMatrix, qdit_index::DimType )

    new_num_qdits = density_matrix.num_qdits-1
    new_dims = zeros(Int8,new_num_qdits)

    skip=0
    for b=1:new_num_qdits
        if b==qdit_index
            skip = 1
            new_dims[b] = density_matrix.dims[b+skip]
            continue
        else
            new_dims[b] = density_matrix.dims[b+skip]
        end
    end

    new_prod_ds = prod_ds(new_dims)

    new_D = prod(new_dims[b] for b=1:new_num_qdits )
    new_rho = zeros(Complex{Float64}, new_D , new_D )


    for I=0:(new_D-1),J=0:(new_D-1)
        for i_a=0:(old_dims[qdit_index]-1)
            old_I = mod(I,old_prod_ds[qdit_index]) + i_qdit_index*old_prod_ds[qdit_index]
            old_J = mod(J,old_prod_ds[qdit_index]) + i_qdit_index*old_prod_ds[qdit_index]

            if a<=new_num_qdits
                old_I = old_I + sum( fld( mod(I,new_prod_ds[b+1] ) , new_prod_ds[b] )*old_prod_ds[b+1] for b=a:new_num_qdits )
                old_J = old_J + sum( fld( mod(J,new_prod_ds[b+1] ) , new_prod_ds[b] )*old_prod_ds[b+1] for b=a:new_num_qdits )
            end

            new_rho[I+1,J+1] = new_rho[I+1,J+1] + density_matrix.rho[old_I+1,old_J+1]
        end
    end

    DensityMatrix(rho=new_rho,dimensions=new_dims,method=from_user)
end


function partial_trace(density_matrix ::DensityMatrix, qdit_indices ::Vector{DimType} )
    qdit_index = sort(qdit_index,rev=True)
    for qdit_index in qdit_indices
        density_matrix = partial_trace(density_matrix,qdit_index)
    end
    density_matrix
end
