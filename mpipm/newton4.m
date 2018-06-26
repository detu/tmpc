N = 20;
global q;
q = 0.0;

x = zeros(4, N);
x(:, 1) = [1; 2; 3; 4];

for i = 1 : N
    [y, J] = f(x(:, i));
%     det(J)
    
    % Newton step
    d = -J \ y;
    
    % Limit Newton step s.t. [mu, t]=x(3:4)>=0 is satisfied
    mu_t = x(3 : 4, i);
%     assert(all(mu_t >= 0));
    
    delta_mu_t = d(3 : 4);
    alpha = max(min(max(-mu_t ./ delta_mu_t, 0), 1));
    
    x(:, i + 1) = x(:, i) + 0.99 * alpha * d;
%     x(:, i + 1) = x(:, i) + d;
end

x

function [y, J] = f(x)
    global q;
    
    mu = x(3);
    t = x(4);
    
    y = [
        1 + 2 * q * x(1) + 2 * mu * x(1);
        1 + 2 * q * x(2) + 2 * mu * x(2);
        2 - x(1)^2 - x(2)^2 - t;
        mu * t
        ];

    J = [
        2 * q + 2 * mu, 0, 2 * x(1), 0;
        0, 2 * q + 2 * mu, 2 * x(2), 0;
        -2 * x(1), -2 * x(2), 0, -1;
        0, 0, t, mu
        ];
end