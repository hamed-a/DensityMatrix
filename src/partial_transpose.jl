using DensityMatrixModule

function partial_transpose(density_matrix::DensityMatrix, qdit_indices )

    new_rho = zeros(Complex{Float64}, density_matrix.total_dimension , density_matrix.total_dimension )

    qdit_indices_a = [1; qdit_indices]

    for I=0:(density_matrix.total_dimension-1),J=0:(density_matrix.total_dimension-1)
        i = [ fld( mod(I,density_matrix.prod_dimensions[b+1] ) , density_matrix.prod_dimensions[b] ) 
            for b=1:density_matrix.num_qdits]
        j = [ fld( mod(J,density_matrix.prod_dimensions[b+1] ) , density_matrix.prod_dimensions[b] ) 
            for b=1:density_matrix.num_qdits]

        i_p = copy(i)
        j_p = copy(j)

        for b in qdit_indices
            i_p[b]=j[b]
            j_p[b]=i[b]
        end

        old_I = sum(i_p[b]*density_matrix.prod_dimensions[b] for b=1:density_matrix.num_qdits)
        old_J = sum(j_p[b]*density_matrix.prod_dimensions[b] for b=1:density_matrix.num_qdits)


        new_rho[I+1,J+1] = density_matrix.rho[old_I+1,old_J+1]
    end

    construct_DensityMatrix(rho=new_rho,qdit_dimensions=density_matrix.qdit_dimensions,method=from_user)
end
