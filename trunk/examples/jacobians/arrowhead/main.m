% This file uses both MATLAB finite differences as well as adigator in order
% to compute derivatives of the arrowhead function. User can change N to
% see the effects of increasing problem size.
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
N = 100;
numeval = 20;
x = rand(N,1);
% ------------------------------ ADiGator ------------------------------- %
gx = adigatorCreateDerivInput([N, 1],'x'); % Create Deriv Input
genout = adigatorGenJacFile('arrowhead',{gx});
S = genout.JacobianStructure;

tic;
for i = 1:numeval
  [Jac,y] = arrowhead_Jac(x);
end
adigatortime = toc/numeval;


% -------------------------- Finite Differences ------------------------- %
TOL = 1e-8;% Can change this to make numjac more/less accurate
tic
for i = 1:numeval
  dfdx = numjac(@arrowhead4numjac,0,x,y,TOL*ones(N,1),[],0);
end
fdtime = toc/numeval;

display(['Average deriv eval time using Finite Differences:           ',num2str(fdtime)]);
display(['Average deriv eval time using ADiGator:                     ',num2str(adigatortime)]);

xs.f = x;
xs.dx = ones(N,1);
