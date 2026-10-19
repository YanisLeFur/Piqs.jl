using QuantumToolbox
export super_tensor

# Equivalent of reshape(matrix, (dimensions..., dimensions...)), followed by
# permutedims with (order..., order .+ N...) and a reshape back to a square matrix,
# done on the indices of the stored entries only.
function permute_dimensions(matrix::SparseMatrixCSC, dimensions, order)
    N = length(dimensions)
    shape = (dimensions..., dimensions...)
    perm = (order..., (order .+ N)...)
    new_shape = map(p -> shape[p], perm)
    return _permute_entries(matrix, CartesianIndices(shape), LinearIndices(new_shape), perm, prod(dimensions))
end

function _permute_entries(matrix, old_indices, new_indices, perm, n)
    rows, cols, vals = findnz(matrix)
    nrows = size(matrix, 1)
    new_linear = map(rows, cols) do i, j
        idx = Tuple(old_indices[i+(j-1)*nrows])
        new_indices[map(p -> idx[p], perm)...]
    end
    return sparse(mod1.(new_linear, n), cld.(new_linear, n), vals, n, n)
end

function compute_order(dimensions)
    shift_vals = 0:2:2*Int(length(dimensions) // 2)-1
    perm_idxs = Tuple(vcat(shift_vals .+ 1, shift_vals .+ 2))
    return perm_idxs
end

"""
    super_tensor(to_tensor)

Tensor product of superoperators.

Given superoperators `L₁, L₂, …, Lₙ` acting on density matrices of the Hilbert
spaces `H₁, H₂, …, Hₙ`, return the superoperator `L` acting on `H₁ ⊗ H₂ ⊗ … ⊗ Hₙ` such that

```math
L\\,[\\rho_1 \\otimes \\rho_2 \\otimes \\dots \\otimes \\rho_n] = L_1[\\rho_1] \\otimes L_2[\\rho_2] \\otimes \\dots \\otimes L_n[\\rho_n].
```

The subsystems may have different dimensions. This is the equivalent of
`super_tensor` in QuTiP.

# Arguments
- `to_tensor`: a tuple (or vector) of `SuperOperator` `QuantumObject`s.

# Returns
- A `SuperOperator` `QuantumObject` whose dimensions are the concatenation of the
  dimensions of the input superoperators.

# Example
```julia
using PIQSolver, QuantumToolbox
a = destroy(2)
b = destroy(3)
L_a = liouvillian(a' * a, [a])
id_b = sprepost(qeye(3), qeye(3))
L = super_tensor((L_a, id_b))   # acts on a 2 ⊗ 3 Hilbert space
```
"""
function super_tensor(to_tensor)
    # the dims of a superoperator are ((row, col), (row, col)) Hilbert-space dims
    hilbert_dims = map(op -> collect(op.dims[1][1]), to_tensor)
    dims = Tuple(reduce(vcat, hilbert_dims))
    # kron and reshape are column-major: the last operator varies fastest, so the
    # (row, column) index pairs of the vectorized density matrix come in reverse order.
    dimensions = tuple((d for hd in reverse(collect(hilbert_dims)) for d in ntuple(_ -> prod(hd), 2))...)
    L_tens = kron(map(op -> op.data, to_tensor)...)
    order = compute_order(dimensions)
    L_tensored = permute_dimensions(L_tens, dimensions, order)
    return Qobj(L_tensored, dims=((dims, dims), (dims, dims)), type=SuperOperator())
end
