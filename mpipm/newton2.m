N = 10;
x = zeros(3, N);
x(:, 1) = -[1; 2; 3];

for i = 1 : N
    % Newton step
    d = -Jf(x(:, i)) \ f(x(:, i));
    
    % Limit Newton step s.t. mu=x(3)>=0 is satisfied
    % ???
    
    x(:, i + 1) = x(:, i) + d;
end

function y = f(x)
    y = [
        1 + 2 * x(3) * x(1);
        1 + 2 * x(3) * x(2);
        x(3) * (2 - x(1)^2 - x(2)^2)
        ];
end

function J = Jf(x)
    J = [
        2 * x(3), 0, 2 * x(1);
        0, 2 * x(3), 2 * x(2);
        -2 * x(3) * x(1), -2 * x(3) * x(2), 2 - x(1)^2 - x(2)^2
        ];
end