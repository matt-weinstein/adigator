function sol = main(n,tspan)
% solves the burgers ode.. ADiGator differentiates the function with the
% loops removed, burgersfun_noloop as it is more efficient for AD
% inputs: n - problem dimension (default of 2.^5)
%         tspan - time span (default is [0 2])
% outputs: sol - solution
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
if nargin == 0
  n = 2.^5;
end
if nargin < 2
  tspan = [0 2];
end

N = n/2;
h = 1/(N+1);
xinit = h*(1:N);
% u(x,0) at grid points
ainit = sin(2*pi*xinit) + 0.5*sin(pi*xinit);

y0 = [ainit xinit];

% Generate ADiGator Jacobian files using adigatorGenJacFile
ay = adigatorCreateDerivInput(size(y0),'y');
output = adigatorGenJacFile('burgersfun_noloop',{1,ay,N},adigatorOptions('overwrite',1));
Jpat = output.JacobianStructure;

% solve ODE
opts = odeset('Mass',@(t,y)burgersmass(t,y,N),'MStateDependence','strong','JPattern',Jpat,...
   'MvPattern',burgersMvPat(N),'RelTol',1e-5,'AbsTol',1e-4,'Jacobian',@(t,y)burgersfun_noloop_Jac(t,y,N));
 display(' ');
sol = ode15s(@(t,y)burgersfun(t,y,N),tspan,y0,opts);

end

function S = burgersMvPat(N)
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.
% Sparsity pattern for the derivative of the Mass matrix times a vector
S = sparse(2*N,2*N);
S(1,2) = 1;
S(1,2+N) = 1;
for i = 2:N-1
   S(i,i-1) = 1;
   S(i,i+1) = 1;
   S(i,i-1+N) = 1;
   S(i,i+1+N) = 1;
end
S(N,N-1) = 1;
S(N,N-1+N) = 1;
end

function out = burgersmass(t,y,N)
%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.
% Mass matrix function. N is provided by the outer function.
a = y(1:N);
x = y(N+1:end);
x0 = 0;
a0 = 0;
xNP1 = 1;
aNP1 = 0;
M1 = speye(N);
M2 = sparse(N,N);
M2(1,1) = - (a(2) - a0)/(x(2) - x0);
for i = 2:N-1
  M2(i,i) = - (a(i+1) - a(i-1))/(x(i+1) - x(i-1));
end
M2(N,N) = - (aNP1 - a(N-1))/(xNP1 - x(N-1));
% MMPDE6
M3 = sparse(N,N);
e = ones(N,1);
M4 = spdiags([e -2*e e],-1:1,N,N);
out = [M1 M2
  M3 M4];
end

