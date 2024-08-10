% Simpson
clear
x = 0;
t = 0;
f = @(s) (s .* sin(s/2 - t) .* cos(x - s));
a = 0;
b = 1;

Iter = 5:1:41;

I = zeros(size(Iter));

for k = Iter
    
    N = k;
    s = linspace(a, b, N);
    vec = f(s);
    h = (b - a)/(N-1);
    I(k-4) = h/3 * (vec(1) + 2 * sum(vec(3:2:end-2)) + 4 * sum(vec(2:2:end)) + vec(end));
    
end
I = abs(I - integral(f, a, b));
loglog(Iter, I )