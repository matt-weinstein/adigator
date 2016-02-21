function [y,varargout] = max(x,varargin)
% CADA overloaded version of function MAX
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if nargout == 2
  error('max only coded for single output')
end
if ADIGATOR.EMPTYFLAG
  y = cadaEmptyEval(x,varargin{:});
  return
end
NUMvod = ADIGATOR.NVAROFDIFF;
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);

% ----------------------------Parse Inputs------------------------------- %
if nargin == 1
  xMrow = x.func.size(1); xNcol = x.func.size(2);
  if xMrow == 1; Dim = 2; else Dim = 1; end
elseif nargin == 2
  y = cadabinaryarraymath(x,varargin{1},0,0,'max');
  return
elseif nargin == 3
  if isa(x,'cada')
    xMrow = x.func.size(1); xNcol = x.func.size(2);
  elseif isnumeric(x)
    [xMrow,xNcol] = size(x);
    xtemp.id = [];
    xtemp.func = struct('name',[],'size',[xMrow xNcol],'zerolocs',[],...
      'value',x);
    if PFLAG
      if xMrow*xNcol == 1
        xtemp.func.name = num2str(x,16);
      else
        xtemp.func.name = cadamatprint(x);
      end
    end
    xtemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    x = xtemp;
  else
    error('??? Cannot sum a non-numeric object.')
  end
  x2  = varargin{1};
  if isa(x2,'cada') || ~isempty(x2)
    error('three input syntax for max: y = max(x,[],dim) -- second input is wrong');
  end
  Dim = varargin{2};
  if isa(Dim,'cada')
    if ~isempty(Dim.func.value)
      ADIGATOR.VARINFO.LASTOCC(Dim.id,1) = ...
        ADIGATOR.VARINFO.COUNT;
      Dim = Dim.func.value;
    else
      error('??? cannot find max across a purely symbolic dimension')
    end
  elseif ~isnumeric(varargin{1})
    error('??? Invalid input argument.')
  end
else
  error('Too many input arguments')
end

if Dim == 1
  yMrow = 1; yNcol = xNcol;
else
  yMrow = xMrow; yNcol = 1;
end
% -----------------------Build Y Function-------------------------------- %
y.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
y.func = struct('name',funcstr,'size',[yMrow yNcol],'zerolocs',[],'value',[]);
if isinf(xMrow) && Dim == 2
  xMrow = 1; xvec = 1;
elseif isinf(xNcol) && Dim == 1
  xNcol = 1; xvec = 2;
elseif isinf(xMrow) || isinf(xNcol)
  error('Cannot find max over vectorized dimension')
else
  xvec = 0;
end


% Function Numerics/Sparsity
if ~isempty(x.func.value)
  % Y is numeric
  y.func.value = max(x.func.value,Dim);
end

% --------------------------Build Derivative----------------------------- %
funprintflag = 0;
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if cadaCheckForDerivs(x)
  % x has derivatives - going to create a dummy variable and call
  % cadamtimesderiv
  TF1 = ['cada',NDstr,'tf1'];
  if DPFLAG
    funprintflag = 1;
    fprintf(fid,[indent,TF1,' = max(',x.func.name,',[],%1.0f);\n'],Dim);
  end
  TF2 = ['cada',NDstr,'tf2'];
  TF3 = ['cada',NDstr,'tf3'];
  TD3 = ['cada',NDstr,'td3'];
  if Dim == 1
    % finding max over first dimension
    A.func.size = [1 xMrow];
    if xNcol == 1 && ~xvec
      % y ~ A*x
      % x is a column vector - make A a row vector of ones and zeros
      if DPFLAG
        fprintf(fid,[indent,TF2,' = (',x.func.name,' == ',TF1,').'';\n']);
        A.func.name = TF2;
      end
    else
      % y ~A*(B.*x)
      % x is a matrix size MxN, make B an MxN matrix of ones and zeros,
      % make A a row vector of ones
      if DPFLAG
        fprintf(fid,[indent,TF2,' = ones(1,%1.0f);\n'],xMrow);
        fprintf(fid,[indent,TF3,' = (',x.func.name,' == repmat(',TF1,',%1.0f,1));\n'],xMrow);
      end
    end
  else
    % finding max over second dimension
    A.func.size = [xNcol 1];
    if xMrow == 1 && ~xvec
      % y ~ x*A
      % x is a row vector, make A a column vector of ones and zeros
      if DPFLAG
        fprintf(fid,[indent,TF2,' = (',x.func.name,' == ',TF1,').'';\n']);
      end
    else
      % x is a matrix size MxN, make B an MxN matrix of ones and zeros,
      % make A a column vector of ones
      if DPFLAG
        fprintf(fid,[indent,TF2,' = ones(%1.0f,1);\n'],xNcol);
        fprintf(fid,[indent,TF3,' = (',x.func.name,' == repmat(',TF1,',1,%1.0f));\n'],xNcol);
      end
    end
  end
  Atemp = ones(A.func.size);
  A.func.name = TF2;
  A.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
  xtemp = ones(xMrow,xNcol);
  if DPFLAG
    fprintf(fid,[indent,funcstr,' = ',TF1,';\n']);
  end
end
for Vcount = 1:NUMvod
  if ~isempty(x.deriv(Vcount).nzlocs)
    derivstr = cadadername(funcstr,Vcount);
    if Dim == 1
      % y = A*X
      if xNcol > 1 || xvec
        % let C = B.*X
        xrows = x.deriv(Vcount).nzlocs(:,1);
        if isequal(xrows,(1:xMrow*xNcol).')
          if xvec
            % Second dim has to be vectorized
            fprintf(fid,[indent,TD3,' = ',TF3,'.''.*',x.deriv(Vcount).name,';\n']);
          else
            fprintf(fid,[indent,TD3,' = ',TF3,'(:).*',x.deriv(Vcount).name,';\n']);
          end
        else
          Dind1 = cadaindprint(xrows);
          if xvec
            % Second dim has to be vectorized
            fprintf(fid,[indent,TD3,' = ',TF3,'(',Dind1,',:).''.*',x.deriv(Vcount).name,';\n']);
          else
            fprintf(fid,[indent,TD3,' = ',TF3,'(',Dind1,').*',x.deriv(Vcount).name,';\n']);
          end
        end
        % Change x deriv name
        x.deriv(Vcount).name = TD3;
      end
      if xvec
        nzlocs = cadamtimesderivvec(A,x,Atemp,xtemp,Vcount,derivstr,DPFLAG,xvec);
      else
        nzlocs = cadamtimesderiv(A,x,Atemp,xtemp,Vcount,derivstr,DPFLAG,'mtimes');
      end
    else
      % y = X*A
      if xMrow > 1 || xvec
        xrows = x.deriv(Vcount).nzlocs(:,1);
        if isequal(xrows,(1:xMrow*xNcol).')
          if xvec
            % First dim has to be vectorized
            fprintf(fid,[indent,TD3,' = ',TF3,'.*',x.deriv(Vcount).name,';\n']);
          else
            fprintf(fid,[indent,TD3,' = ',TF3,'(:).*',x.deriv(Vcount).name,';\n']);
          end
        else
          Dind1 = cadaindprint(xrows);
          if xvec
            % First dim has to be vectorized
            fprintf(fid,[indent,TD3,' = ',TF3,'(:,',Dind1,').*',x.deriv(Vcount).name,';\n']);
          else
            fprintf(fid,[indent,TD3,' = ',TF3,'(',Dind1,').*',x.deriv(Vcount).name,';\n']);
          end
        end
        % Change x deriv name
        x.deriv(Vcount).name = TD3;
      end
      if xvec
        nzlocs = cadamtimesderivvec(x,A,xtemp,Atemp,Vcount,derivstr,DPFLAG,xvec);
      else
        nzlocs = cadamtimesderiv(x,A,xtemp,Atemp,Vcount,derivstr,DPFLAG,'mtimes');
      end
    end
    y.deriv(Vcount).nzlocs = nzlocs;
    y.deriv(Vcount).name = derivstr;
  end
end


% ---------------------Print Out Function-------------------------------- %
if PFLAG && ~funprintflag
  if nargin == 1
    fprintf(fid,[indent,funcstr,' = max(',x.func.name,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = max(',x.func.name,',[],%1.0d);\n'],Dim);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = cada(y);
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;