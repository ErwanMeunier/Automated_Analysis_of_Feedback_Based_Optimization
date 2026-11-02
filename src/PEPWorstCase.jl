module PEPWorstCase

using LinearAlgebra
using ..PEPControl

const PD = PEPControl.PEPDesign

export factor_gram, extract_linear_operator, extract_control_operators,
       solve_discrete_lyapunov, corollary_step_size,
       theorem1_uniform_step_size

"""
    factor_gram(G; rtol=1e-8, atol=1e-10)

Factor a numerical PEP Gram matrix as `G ≈ coordinates' * coordinates` after
discarding eigenvalues below `atol + rtol * λmax`. The returned coordinates
place every PEP basis vector in a common finite-dimensional realization.
"""
function factor_gram(G::AbstractMatrix{<:Real};
                     rtol::Float64 = 1e-8, atol::Float64 = 1e-10)
    rtol >= 0 || error("Gram rank tolerance must be nonnegative")
    atol >= 0 || error("Gram absolute tolerance must be nonnegative")
    decomposition = eigen(Symmetric(Matrix{Float64}(G)))
    eigenvalues = decomposition.values
    lambda_max = max(maximum(eigenvalues; init = 0.0), 0.0)
    threshold = atol + rtol * lambda_max
    keep = findall(>(threshold), eigenvalues)
    order = sort(keep; by = i -> eigenvalues[i], rev = true)
    factor = isempty(order) ? zeros(0, size(G, 1)) :
             Diagonal(sqrt.(eigenvalues[order])) * decomposition.vectors[:, order]'
    approximation = factor' * factor
    relative_error = norm(approximation - G) / max(norm(G), eps(Float64))
    return (; factor, rank = length(order), threshold,
            eigenvalues = sort(eigenvalues; rev = true),
            min_eigenvalue = minimum(eigenvalues; init = 0.0),
            relative_error)
end

function point_coefficients(point, eta, dimension)
    coefficients = zeros(dimension)
    for (index, coefficient) in point.c
        coefficients[index] = PD.evaluate(coefficient, eta)
    end
    return coefficients
end

point_coordinates(point, eta, factor) =
    factor * point_coefficients(point, eta, size(factor, 2))

"""
    extract_linear_operator(operator, eta, factor)

Recover one matrix consistent with every forward pair `(x,Mx)` and adjoint
pair `(u,M'u)` represented in a PEP solution. A minimum-norm subspace
completion selects a deterministic extension on directions not queried by the
PEP, avoiding arbitrary null-space components from an underdetermined solve.
"""
function extract_linear_operator(operator, eta, factor)
    dimension = size(factor, 1)
    dimension > 0 || error("cannot recover an operator from a rank-zero Gram matrix")
    forward = collect(zip(operator.core.xs, operator.core.gs))
    adjoint = collect(zip(operator.adj.xs, operator.adj.gs))
    isempty(forward) && isempty(adjoint) &&
        return zeros(dimension, dimension), 0.0

    X = isempty(forward) ? zeros(dimension, 0) :
        hcat((point_coordinates(pair[1], eta, factor) for pair in forward)...)
    Y = isempty(forward) ? zeros(dimension, 0) :
        hcat((point_coordinates(pair[2], eta, factor) for pair in forward)...)
    U = isempty(adjoint) ? zeros(dimension, 0) :
        hcat((point_coordinates(pair[1], eta, factor) for pair in adjoint)...)
    V = isempty(adjoint) ? zeros(dimension, 0) :
        hcat((point_coordinates(pair[2], eta, factor) for pair in adjoint)...)

    function basis_and_pinv(data)
        isempty(data) && return zeros(dimension, 0), zeros(size(data, 2), dimension)
        decomposition = svd(data; full = false)
        threshold = max(sqrt(eps(Float64)), max(size(data)...) * eps(Float64)) *
                    maximum(decomposition.S; init = 0.0)
        rank = count(>(threshold), decomposition.S)
        rank == 0 && return zeros(dimension, 0), zeros(size(data, 2), dimension)
        basis = decomposition.U[:, 1:rank]
        inverse = decomposition.V[:, 1:rank] *
                  Diagonal(1.0 ./ decomposition.S[1:rank]) * basis'
        return basis, inverse
    end

    Qx, pinv_X = basis_and_pinv(X)
    Qu, pinv_U = basis_and_pinv(U)
    forward_on_basis = Y * pinv_X * Qx
    adjoint_on_basis = V * pinv_U * Qu
    projector_x = Qx * Qx'
    unconstrained = forward_on_basis * Qx' +
                    Qu * adjoint_on_basis' * (I - projector_x)

    # Respect the operator class even when tiny Gram-factor errors make the
    # interpolation system only approximately consistent. Singular-value
    # clipping is the Euclidean projection onto ||M||_2 <= L.
    if iszero(operator.L)
        matrix = zeros(dimension, dimension)
    else
        decomposition = svd(unconstrained)
        matrix = decomposition.U *
                 Diagonal(min.(decomposition.S, operator.L)) *
                 decomposition.Vt
    end

    forward_error = norm(matrix * X - Y)^2
    adjoint_error = norm(matrix' * U - V)^2
    data_scale = norm(Y)^2 + norm(V)^2
    input_scale = operator.L^2 * (norm(X)^2 + norm(U)^2)
    residual = sqrt(forward_error + adjoint_error) /
               max(sqrt(data_scale + input_scale), 1e-10)
    return matrix, residual
end

"""Recover the operators named `A` and `B` from a control PEP model."""
function extract_control_operators(model, eta, factor)
    operators = filter(operator -> operator isa PD.LinearOperator, model.funcs)
    by_name = Dict(operator.core.name => operator for operator in operators)
    haskey(by_name, "A") || error("control PEP contains no linear operator named A")
    haskey(by_name, "B") || error("control PEP contains no linear operator named B")
    A, residual_A = extract_linear_operator(by_name["A"], eta, factor)
    B, residual_B = extract_linear_operator(by_name["B"], eta, factor)
    return (; A, B, residual_A, residual_B,
            bound_A = by_name["A"].L, bound_B = by_name["B"].L)
end

"""
    solve_discrete_lyapunov(A; rtol=1e-11, atol=1e-12, max_iterations=100)

Solve `A'PA - P = -I` by Smith doubling. This is the `ξ=1` Lyapunov
normalization used in the paper's analytical corollary.
"""
function solve_discrete_lyapunov(A::AbstractMatrix{<:Real};
                                 rtol::Float64 = 1e-11,
                                 atol::Float64 = 1e-12,
                                 max_iterations::Int = 100)
    size(A, 1) == size(A, 2) || error("A must be square")
    n = size(A, 1)
    spectral_radius = maximum(abs, eigvals(Matrix{Float64}(A)); init = 0.0)
    spectral_radius < 1 ||
        error("recovered A is not Schur stable (spectral radius=$spectral_radius)")

    power = Matrix{Float64}(A)
    P = Matrix{Float64}(I, n, n)
    converged = false
    iterations = 0
    for iteration in 1:max_iterations
        correction = power' * P * power
        next_P = P + correction
        iterations = iteration
        if norm(correction) <= atol + rtol * max(norm(next_P), 1.0)
            P = next_P
            converged = true
            break
        end
        P = next_P
        power = power * power
    end
    converged || error("discrete Lyapunov doubling did not converge")
    P = Matrix(Symmetric(P))
    residual = norm(A' * P * A - P + I) / max(norm(P), 1.0)
    return (; P, residual, spectral_radius, iterations)
end

"""
    theorem1_uniform_step_size(norm_A, norm_B, L_f)

Return the sharp uniform Proposition-8.3 threshold from Theorem 1,

`2(1-norm_A)^3 / (L_f * norm_B^2 * (3-norm_A))`,

for `0 <= norm_A < 1`, `norm_B > 0`, and `L_f > 0`. The theorem states a
strict condition below this supremal threshold; experiments evaluate the
displayed boundary value itself as the named Theorem-1 comparison point.
"""
function theorem1_uniform_step_size(norm_A::Real, norm_B::Real, L_f::Real)
    a, b, L = Float64(norm_A), Float64(norm_B), Float64(L_f)
    isfinite(a) && 0 <= a < 1 ||
        error("Theorem 1 requires 0 <= norm_A < 1")
    isfinite(b) && b > 0 || error("Theorem 1 requires norm_B > 0")
    isfinite(L) && L > 0 || error("Theorem 1 requires L_f > 0")
    return 2 * (1 - a)^3 / (L * b^2 * (3 - a))
end

"""
    corollary_step_size(A, B, L_f)

Evaluate the paper's discrete-time analytical bound on recovered matrices:
`G=(I-A)\\B`, `L_G=L_f*||G||^2`, `L_tilde=L_f*||G||`,
`zeta=2||A'PG||`, `omega=||G'PG||`, and
`alpha=2/(L_G + L_tilde*(zeta + sqrt(zeta^2 + 4omega)))`.
"""
function corollary_step_size(A::AbstractMatrix{<:Real},
                             B::AbstractMatrix{<:Real}, L_f::Real)
    size(A, 1) == size(A, 2) || error("A must be square")
    size(B, 1) == size(A, 1) || error("A and B have incompatible dimensions")
    isfinite(L_f) && L_f > 0 || error("L_f must be finite and positive")
    lyapunov = solve_discrete_lyapunov(A)
    n = size(A, 1)
    G = (Matrix{Float64}(I, n, n) - A) \ B
    norm_G = opnorm(G)
    L_G = Float64(L_f) * norm_G^2
    L_tilde = Float64(L_f) * norm_G
    zeta = 2 * opnorm(A' * lyapunov.P * G)
    omega = opnorm(G' * lyapunov.P * G)
    denominator = L_G + L_tilde * (zeta + sqrt(zeta^2 + 4omega))
    isfinite(denominator) && denominator > 0 ||
        error("the analytical corollary produced a nonpositive denominator")
    alpha = 2 / denominator
    return (; alpha, G, P = lyapunov.P, L_G, L_tilde, zeta, omega,
            norm_G, lyapunov_residual = lyapunov.residual,
            spectral_radius_A = lyapunov.spectral_radius,
            lyapunov_iterations = lyapunov.iterations)
end

end
