function F = bananaf(x)
% Evaluate the vector function for the system of nonlinear equations
% derived from the general n-dimensional Rosenbrock function.
%   Copyright 1990-2008 The MathWorks, Inc.
% Get the problem size
n = length(x);  
% Evaluate the vector function
odds  = 1:2:n;
evens = 2:2:n;
F = zeros(n,1);
F(odds,1)  = 1-x(odds);
F(evens,1) = 10.*(x(evens)-x(odds).^2); 