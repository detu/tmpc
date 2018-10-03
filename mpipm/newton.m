N = 10;
x = zeros(2, N);
x(:, 1) = [20; 50];

for i = 1 : N
    x(:, i + 1) = x(:, i) - Jf(x(:, i)) \ f(x(:, i));
end

function y = f(x)
    y = [
        x(1) * x(2);
        x(1) + x(2) - 1
        ];
end

function J = Jf(x)
    J = [
        x(2), x(1);
        1, 1
        ];
end