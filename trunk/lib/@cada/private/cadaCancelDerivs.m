function cadaCancelDerivs(varid,numvars)
% this function is called by operations which operate on cada objects which
% may contain derivatives, but the operation always results in an empty
% derivative, 
% ex: a == x(1:3), dont want to print derivs of x(1:3)
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ~ADIGATOR.RUNFLAG
  depvars = find(ADIGATOR.VARINFO.LASTOCC == varid);
  depvars = depvars(depvars > numvars);
  for Vcount = 1:length(depvars)
    depvar = depvars(Vcount);
    if depvar ~= varid
      ADIGATOR.VARINFO.NAMELOCS(depvar,3) = Inf;
      cadaCancelDerivs(depvar,numvars);
    end
  end
end