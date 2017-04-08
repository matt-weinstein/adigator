function [c,S] = adigatorColor(J)
% function [c,S] = adigatorColor(J)
%
% CPR Coloring algorithm - very simplistic.
% You can use coloring if you know your Jacobian sparsity pattern prior to
% calling ADiGator.
%
% J: Jacobian sparsity pattern
%
% c: colors
% S: Seed matrix
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% see also adigatorUncompressJac adigatorCreateDerivInput

n = size(J,2);
c = zeros(1,n);

JJ = J.'*J;
% JJ symmetric matrix, if JJ(i,j) ~= 0, then columns i and j of J share a
% non-zero row value.
x = 1:n;

nc = 0;
for i = 1:n
  % Loop on columns 2 through n
  ci = 0;
  Ji = J(:,i);
  if any(Ji)
    for k = 1:nc
      % Determine if column i can use any of the previously used colors
      xci = x(c==k);
      % xci is collection of columns using color k
      if ~any(JJ(xci,i))
        % if any element of JJ(xci,i) non-zero, then column i shares a 
        % non-zero with one of the columns using color k
        ci = k;
        break
      end
    end
    if ci == 0
      % could not find a suitable color
      ci = nc+1;
      nc = ci;
    end
  end
  c(i) = ci;
end

if nargout == 2
  % Want to build seed matrix as well
  x(c==0) = [];
  c(c==0) = [];
  S = sparse(x,c,1,n,nc);
end