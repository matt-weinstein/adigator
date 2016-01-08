function f = brownf(x)
%BROWNFG Nonlinear minimization problem (function and its gradients).
% Documentation example.

%   Copyright 1990-2008 The MathWorks, Inc.

% Evaluate the function.
n=length(x); y=zeros(n,1)*x(1);
i=1:(n-1);
y(i)=(x(i).^2).^(x(i+1,1).^2+1)+(x(i+1).^2).^(x(i).^2+1);
f=sum(y);