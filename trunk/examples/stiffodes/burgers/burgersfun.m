function out = burgersfun(t,y,N)
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.
% Derivative function.
a = y(1:N);
x = y(N+1:end);
x0 = 0;
a0 = 0;
xNP1 = 1;
aNP1 = 0;
y1 = y(1);
g = zeros(2*N,1).*y1;
for i = 2:N-1
  delx = x(i+1) - x(i-1);
  g(i) = 1e-4*((a(i+1) - a(i))/(x(i+1) - x(i)) - ...
    (a(i) - a(i-1))/(x(i) - x(i-1)))/(0.5*delx) ...
    - 0.5*(a(i+1)^2 - a(i-1)^2)/delx;
end
delx = x(2) - x0;
g(1) = 1e-4*((a(2) - a(1))/(x(2) - x(1)) - (a(1) - a0)/(x(1) - x0))/(0.5*delx) ...
  - 0.5*(a(2)^2 - a0^2)/delx;
delx = xNP1 - x(N-1);
g(N) = 1e-4*((aNP1 - a(N))/(xNP1 - x(N)) - ...
  (a(N) - a(N-1))/(x(N) - x(N-1)))/delx - ...
  0.5*(aNP1^2 - a(N-1)^2)/delx;
% Evaluate monitor function.
M = zeros(N,1).*y1;
for i = 2:N-1
  M(i) = sqrt(1 + ((a(i+1) - a(i-1))/(x(i+1) - x(i-1)))^2);
end
M0 = sqrt(1 + ((a(1) - a0)/(x(1) - x0))^2);
M(1) = sqrt(1 + ((a(2) - a0)/(x(2) - x0))^2);
M(N) = sqrt(1 + ((aNP1 - a(N-1))/(xNP1 - x(N-1)))^2);
MNP1 = sqrt(1 + ((aNP1 - a(N))/(xNP1 - x(N)))^2);
% Spatial smoothing with gamma = 2, p = 2.
SM = zeros(N,1).*y1;
for i = 3:N-2
  SM(i) = sqrt((4*M(i-2)^2 + 6*M(i-1)^2 + 9*M(i)^2 + ...
    6*M(i+1)^2 + 4*M(i+2)^2)/29);
end
SM0 = sqrt((9*M0^2 + 6*M(1)^2 + 4*M(2)^2)/19);
SM(1) = sqrt((6*M0^2 + 9*M(1)^2 + 6*M(2)^2 + 4*M(3)^2)/25);
SM(2) = sqrt((4*M0^2 + 6*M(1)^2 + 9*M(2)^2 + 6*M(3)^2 + 4*M(4)^2)/29);
SM(N-1) = sqrt((4*M(N-3)^2 + 6*M(N-2)^2 + 9*M(N-1)^2 + 6*M(N)^2 + 4*MNP1^2)/29);
SM(N) = sqrt((4*M(N-2)^2 + 6*M(N-1)^2 + 9*M(N)^2 + 6*MNP1^2)/25);
SMNP1 = sqrt((4*M(N-1)^2 + 6*M(N)^2 + 9*MNP1^2)/19);
for i = 2:N-1
  g(i+N) = (SM(i+1) + SM(i))*(x(i+1) - x(i)) - ...
    (SM(i) + SM(i-1))*(x(i) - x(i-1));
end
g(1+N) = (SM(2) + SM(1))*(x(2) - x(1)) - (SM(1) + SM0)*(x(1) - x0);
g(N+N) = (SMNP1 + SM(N))*(xNP1 - x(N)) - (SM(N) + SM(N-1))*(x(N) - x(N-1));
tau = 1e-3;
g(1+N:end) = - g(1+N:end)/(2*tau);
out = g;
end