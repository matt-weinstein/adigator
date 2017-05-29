% This file computes the hessian of the logsumexp problem using
% adigatorGenHesFile
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
clc
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
n = 2^5;

x_x = adigatorCreateDerivInput([n 1],'x');
output = adigatorGenHesFile('logsumexp',{x_x});

x = rand(n,1);

[H,G,f] = logsumexp_Hes(x);
[G2,f2] = logsumexp_Grd(x);
