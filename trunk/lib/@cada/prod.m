function y = prod(x,varargin)
% CADA overloaded version of function PROD
%
% NOTE: If the input has derivatives, this will only work if it is also a
% vector!!
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
if ADIGATOR.EMPTYFLAG
  if nargin == 1
    y = cadaEmptyEval(x);
  else
    y = cadaEmptyEval(x,varargin{1});
  end
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
    error('??? Cannot prod a non-numeric object.')
  end
  Dim = varargin{1};
  if isa(Dim,'cada')
    if ~isempty(Dim.func.value)
      ADIGATOR.VARINFO.LASTOCC(Dim.id,1) = ...
        ADIGATOR.VARINFO.COUNT;
      Dim = Dim.func.value;
    else
      error('??? cannot prod across a purely symbolic dimension')
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
  error('Cannot prod over vectorized dimension')
else
  xvec = 0;
end
if yMrow == 1 && yNcol == 1 && ~xvec
  vectorflag = 1;
else
  vectorflag = 0;
end


% Function Numerics/Sparsity
if ~isempty(x.func.value)
  % Y is numeric
  y.func.value = prod(x.func.value,Dim);
elseif ~isempty(x.func.zerolocs)
  % X is sparse
  xtemp = ones(xMrow,xNcol);
  xtemp(x.func.zerolocs) = 0;
  ytemp = prod(xtemp,Dim);
  y.func.zerolocs = find(~ytemp(:));
else
  xtemp = ones(xMrow,xNcol);
end

% --------------------------Build Derivative----------------------------- %
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if vectorflag && length(x.func.zerolocs) > 1
  % more than one known zero in x, no derivatives in y.
  x.deriv = y.deriv;
end
if cadaCheckForDerivs(x)
  if vectorflag
    % x has derivatives - want to create a dummy variable, P, st. DY = P*DX,
    % then call cadamtimesderiv to perform the projection/multiplication
    % x is a vector, p is vector as well.
    % THREE CASES - will print flow control to the file.
    % if nnz(x) > 1
    %    dy = 0
    % else
    %    if nnz(x) = 1
    %       xtemp = x; xtemp(x == 0) = 1; 
    %       p = zeros(1,n);
    %       p(x==0) = prod(xtemp);
    %    else
    %       p = prod(x)./x;
    %    end
    %    dx = p*dx;
    % end
    nx = xMrow*xNcol;
    x2 = x;
    x2.func.size = [nx 1];
    
    TF2 = ['cada',NDstr,'tf2']; % use for xtemp
    TF3 = ['cada',NDstr,'tf3']; % use for p
    TF4 = ['cada',NDstr,'tf4']; % use for y
    p.func.size = [1 nx];
    p.func.name = TF3;
    p.deriv = y.deriv;
    if ~isempty(x.func.zerolocs)
      % One element is known to be zero - all elements which x~=0 => y=0
      ptemp = zeros(1,nx);
      ptemp(x.func.zerolocs) = 1;
    else
      ptemp = ones(1,nx);
    end
    if DPFLAG
      xnnzStr = ['cada',NDstr,'tnnz1'];
      fprintf(fid,[indent,xnnzStr,' = nnz(',x.func.name,');\n']);
      fprintf(fid,[indent,'cadaconditional1 = ',xnnzStr,' > %1.0d;\n'],nx-2);
      fprintf(fid,[indent,'if cadaconditional1\n']);
      indent = [indent,'    '];
      
      if isempty(x.func.zerolocs)
        % Don't need this statement if we know there is a zero in x.
        fprintf(fid,[indent,'cadaconditional1 = ',xnnzStr,' == %1.0d;\n'],nx-1);
        fprintf(fid,[indent,'if cadaconditional1\n']);
        indent = [indent,'    '];
      end
      
      % Build p vector for nnz(x) == n-1
      Tind = ['cada',NDstr,'tind1'];
      fprintf(fid,[indent,TF2,' = ',x.func.name,';\n']);
      fprintf(fid,[indent,Tind,' = ',TF2,' == 0;\n']);
      fprintf(fid,[indent,TF2,'(',Tind,') = 1;\n']);
      fprintf(fid,[indent,TF3,' = zeros(1,%1.0d);\n'],nx);
      fprintf(fid,[indent,TF3,'(',Tind,') = prod(',TF2,');\n']);
      fprintf(fid,[indent,TF4,' = 0;\n']);
      
      if isempty(x.func.zerolocs)
        fprintf(fid,[indent(5:end),'else\n']);
        % Build p vector for nnz(x) == 0
        fprintf(fid,[indent,TF4,' = prod(',x.func.name,');\n']);
        if xNcol == 1
          fprintf(fid,[indent,TF3,' = (',TF4,'./',x.func.name,').'';\n']);
        else
          fprintf(fid,[indent,TF3,' = ',TF4,'./',x.func.name,';\n']);
        end
        
        indent(1:4) = [];
        fprintf(fid,[indent,'end\n']);
      end
      
      % Print out derivative matrix operation
      ADIGATOR.PRINT.INDENT = indent;
      for Vcount = 1:NUMvod
        if ~isempty(x.deriv(Vcount).nzlocs)
          derivstr = cadadername(funcstr,Vcount);
          nzlocs = cadamtimesderiv(p,x2,ptemp,xtemp(:),Vcount,derivstr,DPFLAG,'mtimes');
          y.deriv(Vcount).nzlocs = nzlocs;
          y.deriv(Vcount).name   = derivstr;
        end
      end
      
      fprintf(fid,[indent(5:end),'else\n']);
      % Derivatives for nnz(x) > 1
      for Vcount = 1:NUMvod
        if ~isempty(y.deriv(Vcount).nzlocs)
          fprintf(fid,[indent,y.deriv(Vcount).name,' = zeros(%1.0d,1);\n'],size(y.deriv(Vcount).nzlocs,1));
        end
      end
      fprintf(fid,[indent,TF4,' = 0;\n']);
      indent(1:4) = [];
      fprintf(fid,[indent,'end\n']);
    else
      for Vcount = 1:NUMvod
        dxinds = x.deriv(Vcount).nzlocs;
        if ~isempty(dxinds)
          if ~isempty(x.func.zerolocs)
            dxinds = dxinds(dxinds(:,1)==x.func.zerolocs,:);
          end
          nv = ADIGATOR.VAROFDIFF(Vcount).usize;
          dx = sparse(dxinds(:,1),dxinds(:,2),ones(size(dxinds,1),1),nx,nv);
          dy = sum(dx,1);
          [yrows,ycols] = find(dy);
          if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
          y.deriv(Vcount).nzlocs = [yrows ycols];
          derivstr = cadadername(funcstr,Vcount);
          y.deriv(Vcount).name = derivstr;
        end
      end
    end

  else
    error('prod not coded for array inputs with derivatives - please use a loop')
  end
end
ADIGATOR.PRINT.INDENT = indent;
% ---------------------Print Out Function-------------------------------- %
if PFLAG
  if cadaCheckForDerivs(x) && DPFLAG
    fprintf(fid,[indent,funcstr,' = ',TF4,';\n']);
  elseif nargin == 1
    fprintf(fid,[indent,funcstr,' = prod(',x.func.name,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = prod(',x.func.name,',%1.0d);\n'],Dim);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = cada(y);
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;