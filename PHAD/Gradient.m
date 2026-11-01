function grad = Gradient(fun, x)
% fun: function handle @(x)

    h = 1e-8;
    n = numel(x);
    grad = zeros(n,1);

    for i = 1:n
        xp = x;
        xm = x;

        xp(i) = xp(i) + h;
        xm(i) = xm(i) - h;

        grad(i) = (fun(xp) - fun(xm)) / (2*h);
    end
end
