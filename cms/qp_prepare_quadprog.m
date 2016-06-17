function problem = qp_prepare_quadprog(qp)
    Nt = length(qp.C);
    [Nx, Nz] = size(qp.C{1});
    Nd = size(qp.D{1}, 1);
    NdT = size(qp.D{end}, 1);

    problem.H = blkdiag(qp.H{:});
    problem.f = vertcat(qp.g{:});
    problem.Aeq = zeros(Nx * Nt, Nz * Nt + Nx);
    problem.Aineq = zeros(Nd * Nt + NdT, Nz * Nt + Nx);

    % x_{k+1} = C_k * z_k + c_k
    % => C_k * z_k - x_{k+1} = -c_k
    for k = 1 : Nt
        problem.Aeq((k - 1) * Nx + (1 : Nx), (k - 1) * Nz + (1 : Nz)) = -qp.C{k};
        problem.Aeq((k - 1) * Nx + (1 : Nx), k * Nz + (1 : Nx)) = eye(Nx);
        problem.Aineq((k - 1) * Nd + (1 : Nd), (k - 1) * Nz + (1 : Nz)) = qp.D{k};
    end
    problem.Aineq(end - NdT + 1 : end, end - Nx + 1 : end) = qp.D{end};
    problem.Aineq = [problem.Aineq; -problem.Aineq];
    problem.bineq = [vertcat(qp.dMax{:}); -vertcat(qp.dMin{:})];

    problem.beq = vertcat(qp.c{:});
    problem.lb = vertcat(qp.zMin{:});
    problem.ub = vertcat(qp.zMax{:});
    
    problem.solver = 'quadprog';
    problem.options = optimoptions('quadprog');
end

