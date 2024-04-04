function residuals = makefwlp_plot(A,b,c,xi,eta,numit,usep)
% makefwlp_plot(A,b,c,xi,eta,numit,usep)
% A,b,c specify SEF LP instance
% xi, eta as in the paper
% numit is number of iterations
% usep is 1 to use 1/(2*sqrt(k)) perturbation terms; 0 for original FWLP

[m,n] = size(A);
x = zeros(n,1);
y = zeros(m,1);
Uvals = zeros(numit,1);
priminfeas = zeros(numit,1);
duinfeas = zeros(numit,1);
dugap = zeros(numit,1);
for k = 1 : numit
    prevx = x;
    prevy = y;
    if usep
        r = xi * simp_project(sqrt(k) * (A' * y - c) / xi);
    else
        [mininfeas,j] = min(c - A' * y);
        r = zeros(n,1);
        if mininfeas < 0
            r(j) = xi;
        end
    end
    x = k * x / (k + 1) + r / (k + 1);
    primres = b - A * x;
    if usep
        primressc = primres * sqrt(k);
        s = zeros(m,1);
        for i = 1 : m
            s(i) = max(-eta, min(eta, primressc(i)));
        end
    else
        s = eta * (primres >= 0) - eta * (primres < 0);
    end
    y = k * y / (k + 1) + s / (k + 1);
    U = -r' * (c - A' * prevy) + s' * (b - A * prevx) + c' * prevx - b' * prevy;
    if usep
        fac = 1 / (2 * sqrt(k));
        U = U - norm(r) * fac - norm(s) * fac;
    end
    Uvals(k) = U;
    priminfeas(k) = norm(b - A * prevx, 1);
    duinfeas(k) = max(0, max(A' * prevy - c));
    dugap(k) = c' * prevx - b' * prevy;
end

plot(1:numit, Uvals, '-k', 1:numit, priminfeas, '-r', 1:numit, duinfeas, '-b', ...
    1:numit, dugap, '-g');

residuals = [Uvals(end), priminfeas(end), duinfeas(end), dugap(end)];

end