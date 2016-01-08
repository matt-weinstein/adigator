function y = subsasgnlogical(x,s,b)
% CADA overloaded version of MATLAB function SUBASGNLOGICAL, to be used in
% conjuction with the ADiGator algorithm.
%
% This should be able to handle:
% y(logical) = x(logical);
% y(logical) = scalar;
% y(logical,logical) = x(logical,logical);
% y(logical,ind) = x(logical) - ind here will only be able to change if the
% first dimension is vectorized
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);
ssize = length(s);
bscalarflag = 0;


% -----------------------Parse B/X input--------------------------------%
if isnumeric(b)
  % b is numeric input
  [bMrow,bNcol] = size(b);
  btemp.id = [];
  btemp.func = struct('name',[],'size',[bMrow,bNcol],'zerolocs',[],...
    'value',b);
  if PFLAG
    if bMrow*bNcol == 1
      btemp.func.name = num2str(b,16);
    else
      cadamatprint(b,['cada',NDstr,'temp1']);
      btemp.func.name = ['cada',NDstr,'temp1'];
    end
  end
  btemp.func.size = [bMrow, bNcol];
  btemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  b = btemp;
  b = class(b,'cada');
else
  bMrow = b.func.size(1);
  bNcol = b.func.size(2);
end
if isnumeric(x)
  % x is numeric input
  [xMrow,xNcol] = size(x);
  xtemp.id = [];
  xtemp.func = struct('name',[],'size',[xMrow,xNcol],'zerolocs',[],...
    'value',x);
  if PFLAG
    if xMrow*xNcol == 1
      xtemp.func.name = num2str(x,16);
    else
      cadamatprint(x,['cada',NDstr,'temp1']);
      xtemp.func.name = ['cada',NDstr,'temp1'];
    end
  end
  xtemp.func.size = [xMrow, xNcol];
  xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  x = xtemp;
  x = class(x,'cada');
else
  xMrow = x.func.size(1);
  xNcol = x.func.size(2);
end

% ----------------------- Parse SUBS inputs --------------------------- %
if length(s.subs) == 2
  [s.subs{1},subs1,logicflag1] = parseIndex(s.subs{1},0);
  [s.subs{2},subs2,logicflag2] = parseIndex(s.subs{2},0);
  if isempty(s.subs{1}) && isinf(xMrow); xMrow = 1; end
  if isempty(s.subs{2}) && isinf(xNcol); xNcol = 1; end
else
  [s.subs{1},subs1,logicflag1] = parseIndex(s.subs{1});
  logicflag2 = 0; subs2 = [];
  if isempty(s.subs{1}) && isinf(xMrow); xMrow = 1; end
  if isempty(s.subs{1}) && isinf(xNcol); xNcol = 1; end
end

% ------------------ Check For Logical Assignment --------------------- %
%logicflag1 => 1st assignment index is logical with unknown values
%logicflag2 => 2nd assignment index is logical with unknown values
if (logicflag1 || logicflag2) && (~bMrow || ~bNcol)
  error(['Cannot use logical assignment with ',...
    'unknown logical values to remove elements of an object'])
end
logicerrflag = 0;
if bMrow*bNcol > 1
  % This could probably be done cleaner, but i think this covers all the
  % cases which are not allowed.
  if logicflag1
    if  ~isfield(b.func,'logicref')
      logicerrflag = 1;
    elseif b.func.logicref(1) ~= subs1.id
      if ~isempty(subs2) && prod(subs2.func.size) > 1
        logicerrflag = 1;
      elseif b.func.size(1) ~= 1 || b.func.logicref(2) ~= subs1.id
        logicerrflag = 1;
      end
    end
  end
  if logicflag2
    if  ~isfield(b.func,'logicref')
      logicerrflag = 1;
    elseif b.func.logicref(2) ~= subs2.id
      if prod(subs1.func.size) > 1
        logicerrflag = 1;
      elseif b.func.size(2) ~= 1 || b.func.logicref(1) ~= subs1.id
        logicerrflag = 1;
      end
    end
  end
end
if logicerrflag
  error(sprintf(['May only use logical assignment with unknown logical ',...
    'values if the object being assigned is either a scalar, or was ',...
    'created by a logical reference, where the logical reference ',...
    'index is the same as the logical assignment index.\n Example:\n',...
    'ind = x < a\nx(ind) = b(ind); is valid.\nx(x<a) = b(x<a); is not.'])); %#ok<SPERR>
end

%---------------------------------------------------------------------%
%                      Build Function Properties                      %
%---------------------------------------------------------------------%
y.id = ADIGATOR.VARINFO.COUNT;
if isinf(xMrow)
  ytemp = 1:xNcol;
  xvec = 1; xnvec = 2; ynvec = 2;
  vecDim = ['size(',x.func.name,',1)'];
  if isinf(bMrow) && isequal(s.subs{1},':')
    bnvec = 2; btemp = 1:bNcol;
  elseif isinf(bNcol) && bMrow==1 && isequal(s.subs{1},':')
    bnvec = 1; btemp = 1;
  elseif bMrow*bNcol==1 && ~cadaCheckForDerivs(b) && isequal(s.subs{1},':')
    bnvec = 0; btemp = 1;
  else
    error('Invalid vectorized subsasgn')
  end
  if length(s.subs) == 2
    ytemp(s.subs{2}) = btemp;
    snvec = 2;
  else
    snvec = 1;
  end
  FMrow = Inf; FNcol = length(ytemp);
elseif isinf(xNcol)
  vecDim = ['size(',x.func.name,',2)'];
  ytemp = (1:xMrow).';
  xvec = 2; xnvec = 1; ynvec = 1;
  if isinf(bNcol) && ...
      (length(s.subs)==2 && isequal(s.subs{2},':')) ||...
      (length(s.subs)==1 && isequal(s.subs{1},':'))
    bnvec = 1; btemp = (1:bMrow).';
  elseif isinf(bMrow) && bNcol==1 && ...
      (length(s.subs)==2 && isequal(s.subs{2},':')) ||...
      (length(s.subs)==1 && isequal(s.subs{1},':'))
    bnvec = 2; btemp = 1;
  elseif bMrow*bNcol==1 && ~cadaCheckForDerivs(b) && ...
      (length(s.subs)==2 && isequal(s.subs{2},':')) ||...
      (length(s.subs)==1 && isequal(s.subs{1},':'))
    bnvec = 0; btemp = 1;
  else
    error('Invalid vectorized subsasgn')
  end
  if length(s.subs) == 2
    ytemp(s.subs{1}) = btemp;
  end
  snvec = 1;
  FMrow = length(ytemp); FNcol = Inf;
elseif prod(x.func.size) == 0 && isinf(bMrow)
  xvec = 1;
  vecDim = ['size(',b.func.name,',1)'];
  ytemp  = [];
  bnvec = 2; btemp = 1:bNcol;
  if length(s.subs) == 2 && isequal(s.subs{1},':')
    ynvec = 2; xnvec = 2;
    FMrow = Inf; FNcol = bNcol;
    snvec = 2;
  elseif length(s.subs) == 2 && isequal(s.subs{2},':') && bNcol == 1
    ynvec = 1; xnvec = 1;
    FMrow = 1; FNcol = Inf;
    snvec = 1;
  else
    error('Invalid vectorized subsasgn')
  end
  ytemp(s.subs{:}) = btemp;
elseif prod(x.func.size) == 0 && isinf(bNcol)
  xvec = 1;
  vecDim = ['size(',b.func.name,',2)'];
  ytemp = [];
  bnvec = 1; btemp = 1:bMrow;
  if length(s.subs) == 2 && isequal(s.subs{1},':') && bMrow == 1
    xnvec = 2; ynvec = 2;
    FMrow = Inf; FNcol = 1;
    snvec = 2;
  elseif length(s.subs) == 2 && isequal(s.subs{2},':')
    xnvec = 1; ynvec = 1;
    FMrow = bMrow; FNcol = Inf;
    snvec = 1;
  else
    error('Invalid vectorized subsasgn')
  end
  ytemp(s.subs{:}) = btemp;
else
  xvec = 0; bnvec =0;
  ytemp = zeros(x.func.size);
  btemp = zeros(b.func.size);
  ytemp(s.subs{:}) = btemp;
  [FMrow, FNcol] = size(ytemp);
end

% Check b scalar case to see if need to REPMAT
if bMrow*bNcol ==1
  btemp = ytemp(s.subs{:});
  [bMrow,bNcol] = size(btemp);
  if bMrow*bNcol > 1
    bscalarflag = 1;
    b.func.size = [bMrow, bNcol];
  elseif xvec
    bscalarflag = 1;
  end
end

if xvec == 1
  ytemp = 1:FNcol;
  if length(s.subs)==2
    ytemp = ytemp(s.subs{2});
  end
elseif xvec == 2
  ytemp = (1:FMrow).';
  if length(s.subs)==2
    ytemp = ytemp(s.subs{1});
  end
else
  ytemp = zeros(FMrow,FNcol);
  ytemp(:) = (1:FMrow*FNcol)';
  ytemp = ytemp(s.subs{:});
end
isubs = ytemp(:);
if xvec == 1
  insubs = 1:FNcol;
elseif xvec==2
  insubs = 1:FMrow;
else
  insubs = 1:FMrow*FNcol;
end
nsubs = length(isubs);
insubs(isubs) = [];
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[FMrow,FNcol],'zerolocs',[],...
  'value',[]);



end

function [numeric, overloaded, logicflag] = parseIndex(index,sflag)
global ADIGATOR
if isa(index,'cada')
  overloaded = index;
  if ~isempty(index.func.value)
    numeric    = index.func.value;
    logicflag  = 0;
  elseif index.func.size(1) == 0 || index.func.size(2) == 0
    numeric    = [];
    logicflag  = 0;
  elseif isfield(index.func,'OnetoN')
    numeric = ':';
    logicflag  = 0;
  elseif isfield(index.func,'logical')
    overloaded = index;
    logicflag  = 1;
    if (isinf(index.func.size(1)) && index.func.size(2) == 1) || ...
        (isinf(index.func.size(2)) && index.func.size(1) == 1)
      if isequal(index.func.value,false)
        numeric = [];
        overloaded.func.name = '[]';
        logicflag = 0;
      else
        numeric = ':';
        logicflag = Inf;
      end
    elseif any(isinf(index.func.size))
      if sflag
        numeric = ':';
        logicflag = Inf;
      else
        error(['A vectorized logical reference index must be of dim ',...
          'N by 1 or 1 by N, where N is vectorized dim']);
      end
    elseif ~isempty(index.func.zerolocs)
      numeric = true(index.func.size);
      numeric(index.func.zerolocs) = false;
    else
      numeric = true(prod(overloaded.func.size),1);
    end
  else
    error('Cannot do strictly symbolic referencing/assignment.')
  end
elseif isnumeric(index)
  logicflag = 0;
  numeric = index;
  overloaded.id = [];
  if ADIGATOR.PRINT.FLAG
    overloaded.func.name = cadaindprint(index);
  end
elseif strcmp(index,':')
  logicflag = 0;
  numeric = index;
  overloaded.id = [];
  overloaded.func.name = ':';
else
  error('Invalid reference index.')
end

end