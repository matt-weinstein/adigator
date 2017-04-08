% This solves the Ginzburg-Landau (2-dimensional) superconductivity problem
% using IPOPT, supplying derivatives via ADiGator and the
% adigatorGenFiles4Ipopt command
% Parameters which can be changed:
% nx - determines problem size
% vorum - is a problem constant which must be an integer
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

% Problem setup
clc
vornum = 8;
nx = 2^4;
ny = nx;
n = 4*nx*ny;

auxdata.nx = nx;
auxdata.ny = ny;
auxdata.tkappa = 5;
auxdata.hx = sqrt(vornum/2)*3/nx;
auxdata.hy = sqrt(vornum/2)*3*sqrt(3)/ny;
auxdata.vornum = vornum;


% Setup structure for adigatorGenFiles4Ipopt
setup.numvar = n;
setup.objective = 'gl2f';
setup.auxdata = auxdata;
setup.order = 2;

% adigatorGenFiles4Ipopt generates everything required by ipopt
tic
funcs = adigatorGenFiles4Ipopt(setup);
gentime = toc;

% Get starting point
z0 = gl2st(auxdata);

% Test callback functions
g0 = feval(funcs.gradient,z0);
h0 = feval(funcs.hessian,z0,1,[]);
Hs = feval(funcs.hessianstructure);

v = rand(n,1);

% % Call ipopt 
if exist('ipopt','file')
  options.ipopt.tol = sqrt(eps);
  options.lb = -Inf*ones(n,1);
  options.ub = Inf*ones(n,1);
  [z, info] = ipopt(z0,funcs,options);
end