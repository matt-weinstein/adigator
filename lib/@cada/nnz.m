function y = nnz(x)
% CADA overloaded version of function NNZ
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x);
  return
end
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;


% -----------------------Build Y Function-------------------------------- %
y.id = ADIGATOR.VARINFO.COUNT;
funcstr = cadafuncname();
y.func = struct('name',funcstr,'size',[1 1],'zerolocs',[],'value',[]);

% Function Numerics/Sparsity
if ~isempty(x.func.value)
  % Y is numeric
  y.func.value = nnz(x.func.value);
end

% --------------------------Build Derivative----------------------------- %
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));

% ---------------------Print Out Function-------------------------------- %
if PFLAG
  fprintf(fid,[indent,funcstr,' = nnz(',x.func.name,');\n']);
end
ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = cada(y);
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;