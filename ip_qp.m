function [X, Lam, Mu, T] = ip_qp(qp)
% Solve a QP using interior point metiod.
% qp is a quadprog-style structure.
    H = qp.H;
    f = qp.f;
    
    lb = qp.lb;
    ub = qp.ub;
    
    nx = size(H, 1);
    assert(isequal(size(H), [nx, nx]));
    assert(isequal(size(f), [nx, 1]));
    assert(isequal(size(lb), [nx, 1]));
    assert(isequal(size(ub), [nx, 1]));
    
    if isfield(qp, 'A')
        A = qp.A;
        b = qp.b;
        nl = size(A, 1);
        
        assert(isequal(size(A), [nl, nx]));
        assert(isequal(size(b), [nl, 1]));
    else
        nl = 0;
        A = zeros(0, nx);
        b = zeros(0, 1);
    end
    
    if isfield(qp, 'Aeq')
        Aeq = qp.Aeq;
        beq = qp.beq;
    
        nm = size(Aeq, 1);
        assert(isequal(size(Aeq), [nm, nx]));
        assert(isequal(size(beq), [nm, 1]));
    else
        nm = 0;
        Aeq = zeros(0, nx);
        beq = zeros(0, 1);
    end
    
    nmbar = nm + 2 * nx;
    
    function [y, J] = F(w)
        x = w(1 : nx);
        lam = w(nx + (1 : nl));
        mu = w(nx + nl + (1 : nmbar));
        t = w(nx + nl + nmbar + (1 : nmbar));
        
        Abar = [A; -speye(nx); speye(nx)];
        bbar = [b; -lb; ub];
        
        y = [
            H * x + f + Aeq.' * lam + Abar.' * mu;
            Aeq * x - beq;
            t + Abar * x - bbar;
            mu .* t
            ];
        
        J = [
            H, Aeq.', Abar.', sparse(nx, nmbar);
            Aeq, sparse(nl, nl + nmbar + nmbar);
            Abar, sparse(nmbar, nl + nmbar), speye(nmbar);
            sparse(nmbar, nx + nl), spdiag(t), spdiag(mu)
            ];
    end

    N = 10;
    w = [zeros(nx, 1); 5 * ones(nl + nmbar + nmbar, 1)];
    
    for i = 1 : N
        [y, J] = F(w);
    %     det(J)

        % Newton step
        d = -J \ y;

        % Limit Newton step s.t. [mu, t]=x(3:4)>=0 is satisfied
        ind_mu_t = nx + (1 : nl + nmbar + nmbar);
        mu_t = w(ind_mu_t);
    %     assert(all(mu_t >= 0));

        delta_mu_t = d(ind_mu_t);
        alpha = max(min(max(-mu_t ./ delta_mu_t, 0), 1));

        w = w + 0.99 * alpha * d;
    %     x(:, i + 1) = x(:, i) + d;
    end

    X = w(1 : nx);
    Lam = w(nx + (1 : nl));
    Mu = w(nx + nl + (1 : nmbar));
    T = w(nx + nl + nmbar + (1 : nmbar));
end

function A = spdiag(v)
    assert(isvector(v));
    n = length(v);
    A = sparse(1 : n, 1 : n, v);
end
