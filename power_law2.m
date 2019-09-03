function S = power_law2(x, b)
alpha = x(1);
beta = x(2);
% gama = x(3);
% S(b>0,1) = beta * (b(b>0)*1000).^-alpha;
S(b>0,1) = beta * (b(b>0)).^-alpha;
S(b == 0,1) = 0;