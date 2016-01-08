function y = arrowsum(x)
N = length(x);

y = zeros(N,1);
y(1) = 2*x(1)^2+sum(x.^2);
y(2:N) = x(1)^2+x(2:N).^2;

y = sum(y);