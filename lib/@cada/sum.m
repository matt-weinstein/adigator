function y = sum(x,varargin)
% CADA overloaded version of function SUM
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% 11/25/14: MJW changed this to no longer call cadamtimesderiv for
% non-vectorized - also, no longer perform sum on columns st dx(:,i) = 0
global ADIGATOR

if ADIGATOR.EMPTYFLAG
  if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ADIGATOR.DERNUMBER == 1
    % If x results from intermediate calculations, size will think the
    % derivatives won't need printed - under the alteration of
    % VARINFO.NAMELOCS
    oldNAMELOCS = ADIGATOR.VARINFO.NAMELOCS;
    sumsize = size(x,1);
    ADIGATOR.VARINFO.NAMELOCS = oldNAMELOCS;
  end
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
    error('??? Cannot sum a non-numeric object.')
  end
  Dim = varargin{1};
  if isa(Dim,'cada')
    if ~isempty(Dim.func.value)
      ADIGATOR.VARINFO.LASTOCC(Dim.id,1) = ...
        ADIGATOR.VARINFO.COUNT;
      Dim = Dim.func.value;
    else
      error('??? cannot sum across a purely symbolic dimension')
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
if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ADIGATOR.DERNUMBER == 1
  sumsize = size(x,Dim);
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
  error('Cannot sum over vectorized dimension')
else
  xvec = 0;
end


% Function Numerics/Sparsity
if ~isempty(x.func.value)
  % Y is numeric
  y.func.value = sum(x.func.value,Dim);
elseif ~isempty(x.func.zerolocs)
  % X is sparse
  xtemp = true(xMrow,xNcol);
  xtemp(x.func.zerolocs) = false;
  ytemp = sum(xtemp,Dim);
  y.func.zerolocs = find(~ytemp(:));
end
% --------------------------Build Derivative----------------------------- %
y.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
if xvec
  % VECTORIZED CALLS CADAMTIMESDERIVVEC
  if cadaCheckForDerivs(x)
    % x has derivatives - going to create a dummy variable and call
    % cadamtimesderiv
    TF2 = ['cada',NDstr,'tf2'];
    if Dim == 1
      % A is 1 by xMrow
      A.func.size = [1 xMrow];
      % - sum(x,1) is same as A*X
    else
      % A is xNcol by 1
      A.func.size = [xNcol 1];
      % - sum(x,2) is same as X*A
    end
    A.func.name = TF2;
    if DPFLAG
      fprintf(fid,[indent,TF2,' = ones(%1.0f,%1.0f);\n'],...
        A.func.size(1),A.func.size(2));
    end
    Atemp = ones(A.func.size);
    A.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    xtemp = ones(xMrow,xNcol);
  end
  for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).nzlocs)
      derivstr = cadadername(funcstr,Vcount);
      if DPFLAG && ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ...
          ADIGATOR.DERNUMBER == 1 && isempty(sumsize.func.value)
        % size(x,Dim) changes - make sure that we don't take zero
        % derivatives.
        TD2 = ['cada',NDstr,'td2'];
        Istr = cadaindprint(x.deriv(Vcount).nzlocs(:,1));
        fprintf(fid,[indent,TD2,' = ',x.deriv(Vcount).name,';\n']);
        fprintf(fid,[indent,TD2,'(:,',Istr,'>',sumsize.func.name,') = 0;\n']);
        x.deriv(Vcount).name = TD2;
      end
      if Dim == 1
        % y = A*X
        nzlocs = cadamtimesderivvec(A,x,Atemp,xtemp,Vcount,derivstr,DPFLAG,xvec);
      else
        % y = X*A
        nzlocs = cadamtimesderivvec(x,A,xtemp,Atemp,Vcount,derivstr,DPFLAG,xvec);
      end
      y.deriv(Vcount).nzlocs = nzlocs;
      y.deriv(Vcount).name = derivstr;
    end
  end
else
  % NONVECTORIZED
  for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).nzlocs)
      derivstr = cadadername(funcstr,Vcount);
      xrows = x.deriv(Vcount).nzlocs(:,1);
      xcols = x.deriv(Vcount).nzlocs(:,2);
      nzx = length(xrows);
      nv = ADIGATOR.VAROFDIFF(Vcount).usize;
      
      if (Dim == 1 && xMrow == 1) || (Dim == 2 && xNcol == 1)
        y.deriv(Vcount).nzlocs = x.deriv(Vcount).nzlocs;
        nzy = nzx;
      elseif (Dim == 1 && xNcol == 1) || (Dim == 2 && xMrow == 1)
        % summing a vector
        dx = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
        dy = sum(dx,1);
        [yrows,ycols] = find(dy);
        if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
        xdim = xMrow*xNcol;
        nzy = length(yrows);
        y.deriv(Vcount).nzlocs = [yrows,ycols];
      elseif Dim == 1
        % x is a matrix, sum over first dimension
        dx = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
        dxR = reshape(dx,xMrow,xNcol*nv);
        dyR  = sum(dxR,1);
        dy   = reshape(dyR,xNcol,nv);
        [yrows,ycols] = find(dy);
        if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
        y.deriv(Vcount).nzlocs = [yrows,ycols];
        if DPFLAG
          % If printing, want to change ycols, dx, xrows, xcols, xdim, nv
          % to make it act like the vector case
          [~,ycols] = find(dyR);
          dx = dxR;
          [xrows,xcols] = find(dxR);
          nzy = length(yrows);
          nv = nv*xNcol;
          xdim = xMrow;
        end
      else
        % x is matrix, sum over second dim - need to change dims up
        xref = zeros(xMrow,xNcol); xref(:) = 1:xMrow*xNcol;
        xref = xref.';
        dx   = sparse(xrows,xcols,1:nzx,xMrow*xNcol,nv);
        dxt  = dx(xref(:),:);
        dxtR = reshape(dxt,xNcol,xMrow*nv);
        dyR  = sum(dxtR,1);
        dy   = reshape(dyR,xMrow,nv);
        [yrows,ycols] = find(dy);
        if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
        y.deriv(Vcount).nzlocs = [yrows,ycols];
        if DPFLAG
          % Change ycols, dx, xrows, xcols, xdim, nv
          [~,ycols] = find(dyR);
          dx = dxtR;
          [xrows,xcols] = find(dxtR);
          nzy = length(yrows);
          nv  = nv*xMrow;
          xdim = xNcol;
        end
      end
      y.deriv(Vcount).name = derivstr;
      if DPFLAG
        LoopSizeChange = 0;
        if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ...
            ADIGATOR.DERNUMBER == 1 && isempty(sumsize.func.value)
          % size(x,Dim) changes - make sure that we don't take zero
          % derivatives.
          TD2 = ['cada',NDstr,'td2'];
          Istr = cadaindprint(xrows);
          fprintf(fid,[indent,TD2,' = ',x.deriv(Vcount).name,';\n']);
          fprintf(fid,[indent,TD2,'(',Istr,'>',sumsize.func.name,') = 0;\n']);
          x.deriv(Vcount).name = TD2;
          LoopSizeChange = 1;
        end
        if nzx == nzy
          % This is a special case meaning that there is only a single
          % entry in each column of dx, thus it is not necessary to to
          % even take a sum. Due to the ordering, it is the case that
          % nonzeros(dx) = nonzeros(dy).
          fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
        else
          % Need to Project dx into a matrix
          if nzy < nv
            % dy not full, do not need to project dx into full nx by nv
            % matrix - ycols gives us the non-zero columns of x.
            [xrows,xcols] = find(dx(:,ycols));
          end
          TD1 = ['cada',NDstr,'td1'];
          if xdim*nzy >= 250 && nzx/(xdim*nzy) <= (3/4)
            % Project sparse
            if ~LoopSizeChange; Istr = cadaindprint(xrows); end
            Jstr = cadaindprint(xcols);
            fprintf(fid,[indent,TD1,' = sum(sparse(',Istr,',',Jstr,',',...
              x.deriv(Vcount).name,',%1.0d,%1.0d),1);\n'],xdim,nzy);
            fprintf(fid,[indent,derivstr,' = full(',TD1,'(:));\n']);
          else
            % Project full
            K = (xcols-1)*xdim + xrows;
            Kstr = cadaindprint(K);
            fprintf(fid,[indent,TD1,' = zeros(%1.0d,%1.0d);\n'],xdim,nzy);
            fprintf(fid,[indent,TD1,'(',Kstr,') = ',x.deriv(Vcount).name,';\n']);
            fprintf(fid,[indent,TD1,' = sum(',TD1,',1);\n']);
            fprintf(fid,[indent,derivstr,' = ',TD1,'(:);\n']);
          end
        end
      end
    end
  end
end

% ---------------------Print Out Function-------------------------------- %
if PFLAG
  if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && ...
      ADIGATOR.DERNUMBER == 1 && isempty(sumsize.func.value)
    % Sum dimension changes - make sure not summing values which
    % arent supposed to be there
    TF1 = 'cada1tempf1';
    fprintf(fid,[indent,TF1,' = ',x.func.name,';\n']);
    if Dim == 1
      FunInd = cadaindprint(1:xMrow);
      fprintf(fid,[indent,TF1,'(',FunInd,'>',sumsize.func.name,',:) = 0;\n']);
    else
      FunInd = cadaindprint(1:xNcol);
      fprintf(fid,[indent,TF1,'(:,',FunInd,'>',sumsize.func.name,') = 0;\n']);
    end
    x.func.name = TF1;
  end
  if nargin == 1
    fprintf(fid,[indent,funcstr,' = sum(',x.func.name,');\n']);
  else
    fprintf(fid,[indent,funcstr,' = sum(',x.func.name,',%1.0d);\n'],Dim);
  end
end

ADIGATOR.VARINFO.LASTOCC([x.id y.id],1) = ADIGATOR.VARINFO.COUNT;
y = cada(y);
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;