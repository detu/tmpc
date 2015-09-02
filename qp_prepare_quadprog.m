function problem = qp_prepare_quadprog(qp)
    Nt = length(qp.C);
    [Nx, Nz] = size(qp.C{1});

    problem.H = blkdiag(qp.H{:});
    problem.g = vertcat(qp.g{:});
    problem.Aeq = zeros(Nx * Nt, Nz * Nt + Nx);

    % x_{k+1} = C_k * z_k + c_k
    % => C_k * z_k - x_{k+1} = -c_k
    for k = 1 : Nt
        problem.Aeq((k - 1) * Nx + (1 : Nx), (k - 1) * Nz + (1 : Nz)) = -qp.C{k};
        problem.Aeq((k - 1) * Nx + (1 : Nx), k * Nz + (1 : Nx)) = eye(Nx);
    end

    problem.beq = vertcat(qp.c{:});
    problem.lb = vertcat(qp.zMin{:});
    problem.ub = vertcat(qp.zMax{:});
    
    problem.solver = 'quadprog';
    problem.options = optimoptions('quadprog');
end

