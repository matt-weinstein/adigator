function m = adigatorFindMatchingParen(str,n)
% Finds the parenthesis pair given an open or close location. Modified from
% MATLAB Cody problem to identifier type of bracket e.g. (, [, {
% Inputs:
%   str - string to parse
%     n - location of open/close
% Copyright 2011-214 Matthew J. Weinstein
% Distributed under the GNU General Public License version 3.0

if strcmp(str(n),'{') || strcmp(str(n),'}')
  os = '{'; cs = '}';
elseif strcmp(str(n),'(') || strcmp(str(n),')')
  os = '('; cs = ')';
elseif strcmp(str(n),'[') || strcmp(str(n),']')
  os = '['; cs = ']';
else
  error('String location not a bracket');
end

startLocs = strfind(str,os);
endLocs   = strfind(str,cs);


if any(startLocs == n)
  forward = true;
else
  % If given closing location, just flip problem
  forward = false;
  ll   = 1:length(str);
  str2 = fliplr(str);
  ll2  = fliplr(ll);
  startLocs = strfind(str2,cs);
  endLocs   = strfind(str2,os);
  n    = find(ll2 == n);
end

endLocs = endLocs(endLocs > n);
startLocs = startLocs(startLocs > n);
m     = 0;
for k = 1:length(endLocs)
  nb = nnz(startLocs < endLocs(k));
  if nb == k-1
    m = endLocs(k);
    break;
  end
end

if forward == false
  m = ll2(m);
end
end