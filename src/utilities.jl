using QuantumToolbox
export super_tensor

function permute_dimensions(matrix, dimensions, order)
    reshape_factor = (dimensions..., dimensions...)
    tensor_data = _reshape_sparse(matrix, reshape_factor)
    N = length(dimensions)
    perm = (order..., order .+ N...)
    tensor_data = permutedims(tensor_data, perm)
    new_shape = (prod(dimensions), prod(dimensions))
    return _to_csc(_reshape_sparse(tensor_data, new_shape))
end

function compute_order(dimensions)
    shift_vals = 0:2:2*Int(length(dimensions) // 2)-1
    perm_idxs = Tuple(vcat(shift_vals .+ 1, shift_vals .+ 2))
    return perm_idxs
end

function super_tensor(to_tensor)
    dims = map(op -> op.dims[1], to_tensor)
    dimensions = tuple((d for op_dim in dims for d in (op_dim, op_dim))...)
    L_tens = kron(map(op -> op.data, to_tensor)...)
    order = compute_order(dimensions)
    L_tensored = permute_dimensions(L_tens, dimensions, order)
    return Qobj(L_tensored, dims=dims, type=SuperOperator)
end