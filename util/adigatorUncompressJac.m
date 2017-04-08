function J = adigatorUncompressJac(Jpat,c,JSnz)
% function J = adigatorUncompressJac(Jpat,c,JSnz)
%
% Uncompress Jacobian
%
% Jpat: Jacobian Sparsity Pattern
% c:    colors
% JSnz: Non-zeros of the compressed Jacobian
%
% J:    Jacobian
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% see also adigatorColor adigatorCreateDerivInput

[m,n] = size(Jpat);
[i,j] = find(Jpat);
order = nonzeros(sparse(i,c(j),1:length(i),m,n));
J = sparse(i,j,JSnz(order),m,n);