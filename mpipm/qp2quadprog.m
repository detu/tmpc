function problem = qp2quadprog(qp)
    n_nodes = length(qp.nodes);
    
    for k = n_nodes : -1 : 1
        qp1(k).H = blkdiag([qp(k).Q, qp(k).S; qp(k).S.', qp(k).R], qp(k).Zl, qp(k).Zu);
        qp1(k).g = [qp(k).q; qp(k).r; qp(k).zl; qp(k).zu];
        
        % state lower bound; input lower bound; lower slack lower bound; upper slack lower bound.
        qp1(k).lbz = [qp(k).lbx; qp(k).lbu; zeros(size(qp(k).zl)); zeros(size(qp(k).zu))];
        
        % state upper bound; input upper bound; lower slack upper bound; upper slack upper bound.
        qp1(k).ubz = [qp(k).ubx; qp(k).ubu; inf(size(qp(k).zl)); inf(size(qp(k).zu))];
        
        qp1(k).C = [qp(k).A, qp(k).B];
        qp1(k).c = qp(k).b;
        qp1(k).D = [qp(k).C, qp(k).D];
        qp1(k).lbd = qp(k).lbd;
        qp1(k).ubd = qp(k).ubd;
    end
    
    qp = qp1;
    
    nz = cellfun(@length, {qp.g});
    nx_next = cellfun(@length, {qp.c});
    nd = cellfun(@length, {qp.lbd});
    
    problem.H = blkdiag(qp.H);
    problem.f = vertcat(qp.g);
    problem.Aeq = zeros(sum(nx_next(1 : end - 1)), sum(nz));
    problem.Aineq = zeros(sum(nd), sum(nz));

    % x_{k+1} = C_k * z_k + c_k
    % => C_k * z_k - x_{k+1} = -c_k
    iz = 0;
    ix_next = 0;
    id = 0;
    for k = 1 : n_nodes
        if k < n_nodes
            problem.Aeq(ix_next + (1 : nx_next(k)), iz + (1 : nz(k) + nx_next(k))) = [-qp(k).C, eye(nx_next(k))];
        end
        
        problem.Aineq(id + (1 : nd(k)), iz + (1 : nz(k))) = qp(k).D;
        
        iz = iz + nz(k);
        ix_next = ix_next + nx_next(k);
        id = id + nd(k);
    end
    problem.Aineq = [problem.Aineq; -problem.Aineq];
    problem.bineq = [vertcat(qp.ubd); -vertcat(qp.lbd)];

    problem.beq = vertcat(qp(1 : end - 1).c);
    problem.lb = vertcat(qp.lbz);
    problem.ub = vertcat(qp.ubz);
    
    problem.solver = 'quadprog';
    problem.options = optimoptions('quadprog');
end

