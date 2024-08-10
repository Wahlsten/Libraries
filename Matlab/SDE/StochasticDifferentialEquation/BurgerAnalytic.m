clear
clc
epsilon = .1;
xv = -5:.1:5;
tv = [3];
l = 0;

for t = tv
    
    l = l + 1;
    k = 0;
    
    for x = xv
        
        k = k + 1;
        g       = @(y) (exp(-10*y.^2));
        theta0  = @(y) arrayfun(@(yk) exp(- 1/(2*epsilon) * integral(g, 0, yk)), y);
        f1      = @(y) ((x - y)./t .* exp(-abs(x - y).^2 ./ (4 * epsilon * t)) .* theta0(y));
        f2      = @(y) (exp(-abs(x - y).^2 ./ (4 * epsilon * t)) .* theta0(y));
        
        funk(l, k) = integral(f1, -10 , 10)/integral(f2, -10 , 10);
        
    end
    
end
hold on
plot(xv, funk)
