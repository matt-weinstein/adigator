function nzlocs = cadamtimesderiv(x,y,xtemp,ytemp,Vcount,derivstr,DPFLAG,caller)
% This function does the Derivative calculations for Matrix Multiplication,
% where it may be called from mtimes, mldivide, or mrdivide. We note that
% the multiplication a*db for matrices is relatively simple given the
% unrolled form in which we store derivatives, thus we use the following
% derivative rules. Allowing that dx' is the derivative of the transpose of
% x, x'
%
% If called from mtimes: 
%   z  = x*y
%   dz = dx*y + x*dy = (y'*dx')' + x*dy
%
% If called from mldivide: 
%   z  = x\y and this will only be called with y having derivatives, thus
%   dz = x\dy
%   see mldivide for more details
%
% If called from mrdivide:
%   z  = x/y and this will only be called with x having derivatives, thus
%   dz = dx/y = (y'\dx')'
%   see mrdivide for more details
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% 11/25/14: MJW made it such that we no longer operate on columns of
% derivative matrices which are known to be all zero.

global ADIGATOR
if DPFLAG
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  NDstr  = sprintf('%1.0f',ADIGATOR.DERNUMBER);
end

xMrow = x.func.size(1); xNcol = x.func.size(2);
yMrow = y.func.size(1); yNcol = y.func.size(2);
switch caller
  case 'mtimes'
    FMrow = xMrow;  FNcol = yNcol;
  case 'mldivide'
    FMrow = xNcol;  FNcol = yNcol;
  case 'mrdivide'
    FMrow = xMrow;  FNcol = yMrow;
  otherwise
   error('??? invalid caller specified') 
end


dxind = x.deriv(Vcount).nzlocs; dyind = y.deriv(Vcount).nzlocs;
nzx = size(dxind,1); nzy = size(dyind,1);
nv = ADIGATOR.VAROFDIFF(Vcount).usize;
if ~isempty(dxind) && ~isempty(dyind)
  % ---------------------X and Y have Derivatives---------------------- %
  % derivatives of Z from Y can be found by
  % X*reshape(DY,yMrow,yNcol*nv) where DY is the projection of DY using
  % the subs index in dyind.
  % derivatives of Z from X can be found using a similiar approach,
  % though we say Z.' = Y.'*X.' and find the derivatives of Z.'
  if ~strcmp(caller,'mtimes')
    error('?? caller can only be mtimes if both variables have derivatives')
  end
  
  % Get DZ from Y first.
  dy   = sparse(dyind(:,1),dyind(:,2),1:nzy,yMrow*yNcol,nv); % DY
  dyR  = reshape(dy,yMrow,yNcol*nv); % DY reshaped
  dzyR = xtemp*dyR; % DZDY reshaped
  dzy  = reshape(dzyR,FMrow*FNcol,nv); % DZDY
  [dzy_rows, dzy_cols] = find(dzy);
  if size(dzy_rows,2) > 1;dzy_rows = dzy_rows.';dzy_cols = dzy_cols.';end
  
  % Get DZ from X.
  xTranMap    = zeros(xMrow,xNcol);
  xTranMap(:) = 1:xMrow*xNcol;
  xTranMap    = xTranMap.';
  dx          = sparse(dxind(:,1),dxind(:,2),1:nzx,xMrow*xNcol,nv);
  dxTran      = dx(xTranMap(:),:); % DX.'
  dxTranR     = reshape(dxTran,xNcol,xMrow*nv); % DX' reshaped
  dzxTranR    = ytemp.'*dxTranR;   % DZDX' reshaped
  dzxTran     = reshape(dzxTranR,FMrow*FNcol,nv); % DZDX'
  zTranMap    = zeros(FNcol,FMrow); zTranMap(:) = 1:FMrow*FNcol;
  zTranMap    = zTranMap.';
  dzx         = dzxTran(zTranMap(:),:); % DZDX
  [dzx_rows, dzx_cols] = find(dzx);
  if size(dzx_rows,2) > 1;dzx_rows = dzx_rows.';dzx_cols = dzx_cols.';end
  
  if ~isempty(dzy_rows) || ~isempty(dzx_rows)
    % Can Now build DZ
    dz = sparse([dzy_rows;dzx_rows],[dzy_cols;dzx_cols],...
      [-ones(size(dzy_rows));2*ones(size(dzx_rows))],FMrow*FNcol,nv);
    % DZ with val < 2 corresponds to Y, DZ with val > 0 corresponds with
    % X
    [zrows,zcols] = find(dz);
    if size(zrows,2) > 1; zrows = zrows.'; zcols = zcols.'; end
    nzlocs = [zrows,zcols];
    
    if DPFLAG
      % ----------------Derivative Printing---------------------------- %
      nz = length(zrows);
      % Print out the Calculations for DZX first. -
      TD1 = ['cada',NDstr,'td1'];
      if ~isempty(dzx_rows)
        % We need to project the derivatives of DX into a matrix of the
        % form of reshape(dxTran,xNcol,xMrow*nv) - but remove zeros columns
        TD2 = ['cada',NDstr,'td2'];
        
        % Determine how many columns of DX' reshaped actually have non-zero
        % derivatives
        dxTranRsum = sum(dxTranR,1);
        dxTranRsumlocs = find(dxTranRsum);
        coldim = length(dxTranRsumlocs);
        if coldim < xMrow*nv
          [dxTRrows, dxTRcols, dxTRlocs] = find(dxTranR(:,dxTranRsumlocs));
        else
          [dxTRrows, dxTRcols, dxTRlocs] = find(dxTranR);
        end
        [~,dxTRlocs] = sortrows(dxTRlocs(:));
        dxTRrows = dxTRrows(dxTRlocs);
        dxTRcols = dxTRcols(dxTRlocs);
        
        % Project dx nonzeros into matrix size xN x coldim
        if xNcol*coldim > 250 && nzx < (3/4)*xNcol*coldim
          SPXflag = 1;
          % Do the projection sparsely - dxTran is larger than 250 and
          % has at least 25% zeros.
          TDind1 = cadaindprint(dxTRrows);
          TDind2 = cadaindprint(dxTRcols);
          % Print out the Projection
          fprintf(fid,[indent,TD2,' = sparse(',TDind1,',',TDind2,',',...
            x.deriv(Vcount).name,',%1.0d,%1.0d);\n'],xNcol,coldim);
        else
          SPXflag = 0;
          % Just project into zeros using a linear index - linear index
          % of dxTran and dxTran-reshaped are the same
          TDind1 = cadaindprint(sub2ind([xNcol,coldim],dxTRrows,dxTRcols));
          fprintf(fid,[indent,TD2,' = zeros(%1.0d,%1.0d);\n'],...
            xNcol,coldim);
          fprintf(fid,[indent,TD2,'(',TDind1,') = ',...
            x.deriv(Vcount).name,';\n']);
        end
        
        % Can now print out Y.'*DX.' - where DX.' is in file as TD2.
        fprintf(fid,[indent,TD2,' = ',y.func.name,'.''*',TD2,';\n']);
        
        % Need to get the mapping from dzxTran to dzx.
        if coldim < xMrow*nv
          [dzxTRrows,dzxTRcols] = find(dzxTranR(:,dxTranRsumlocs));
          dzxTranR = sparse(dzxTRrows,dxTranRsumlocs(dzxTRcols),...
            sub2ind([FNcol,coldim],dzxTRrows,dzxTRcols),FNcol,FMrow*nv);
          dzxTran = reshape(dzxTranR,FMrow*FNcol,nv);
        else
          [dzxTrows,dzxTcols] = find(dzxTran);
          dzxTran  = sparse(dzxTrows,dzxTcols,...
            sub2ind([FMrow*FNcol,nv],dzxTrows,dzxTcols),FMrow*FNcol,nv);
        end
        dzx      = dzxTran(zTranMap(:),:);
        dzxTinds = nonzeros(dzx);
        Dind1    = cadaindprint(dzxTinds(:));
        % dzxTinds are now the Linear Reference index off of dzxTran (in
        % the file) to the nonzeros in dzx.

        % need to find which parts of DZ correspond to DZX
        dzxInds = 1:nz;
        dzxInds = dzxInds(nonzeros(dz)>0);
        if length(dzxInds) == nz
          if SPXflag
            fprintf(fid,[indent,TD1,' = full(',TD2,'(',Dind1,'));\n']);
          else
            fprintf(fid,[indent,TD1,' = ',TD2,'(',Dind1,');\n']);
          end
          fprintf(fid,[indent,TD1,' = ',TD1,'(:);\n']);
        else
          Dind2 = cadaindprint(dzxInds);
          fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],nz);
          fprintf(fid,[indent,TD1,'(',Dind2,') = ',TD2,'(',Dind1,');\n']);
        end
      end
      
      if ~isempty(dzy_rows)
        % Print out the Calculations for DZY
        TD2 = ['cada',NDstr,'td2'];
        % Need to project the vector of DY(in the file) into the matrix
        % of the form reshape(DY,yMrow,yNcol*nv) - except remove any known
        % zero columns
        
        % Determine how many columns of DY reshaped actually have non-zero
        % derivatives
        dyRsum = sum(dyR,1);
        dyRsumlocs = find(dyRsum);
        coldim = length(dyRsumlocs);
        if coldim < yNcol*nv
          % dyR has some columns which are all zero - do not need to sum
          % these columns
          [dyRrows, dyRcols] = find(dyR(:,dyRsumlocs));
        else
          [dyRrows, dyRcols] = find(dyR);
        end
        
        % Project dy nonzeros into matrix size yM by coldim
        if yMrow*coldim > 250 && nzy < (3/4)*yMrow*coldim
          SPYflag = 1;
          % Do the projection sparsely - dy is larger than 250 and has at
          % least %25 zeros.
          TDind1 = cadaindprint(dyRrows);
          TDind2 = cadaindprint(dyRcols);
          % Print out the Projection
          fprintf(fid,[indent,TD2,' = sparse(',TDind1,',',TDind2,',',...
            y.deriv(Vcount).name,',%1.0d,%1.0d);\n'],yMrow,coldim);
        else
          SPYflag = 0;
          % Project into zeros using a linear index.
          TDind1 = cadaindprint(sub2ind([yMrow,coldim],dyRrows,dyRcols));
          % Print out the Projection
          fprintf(fid,[indent,TD2,' = zeros(%1.0d,%1.0d);\n'],...
            yMrow,coldim);
          fprintf(fid,[indent,TD2,'(',TDind1,') = ',...
            y.deriv(Vcount).name,';\n']);
        end
        
        % Can now print out X*DY - where DY is TD2
        fprintf(fid,[indent,TD2,' = ',x.func.name,'*',TD2,';\n']);
        fprintf(fid,[indent,TD2,' = ',TD2,'(:);\n']);
        % TD2 in the file now corresponds with DZY (except reshaped), so
        % we need to know what are the nonzero indices of DZY which we
        % need to reference off of TD2 (even though it is reshaped,
        % linear indexing will still be the same)
        
        if coldim < yNcol*nv
          dyInds = find(dzyR(:,dyRsumlocs));
        else
          dyInds = sub2ind([FMrow*FNcol,nv],dzy_rows,dzy_cols);
        end

        Dind1  = cadaindprint(dyInds(:));
        % dyInds are now the linear reference off of TD2 which give the
        % nonzeros in DZY - just need to figure out the mapping of DZY
        % into DZ now.
        dzyInds = 1:nz;
        dzyInds = dzyInds(nonzeros(dz) < 2);

        if length(dzyInds) == nz
          if SPYflag
            if ~isempty(dzx_rows)
              fprintf(fid,[indent,TD1,' = ',TD1,' + full(',TD2,'(',Dind1,'));\n']);
            else
              fprintf(fid,[indent,TD1,' = full(',TD2,'(',Dind1,'));\n']);
            end
          else
            if ~isempty(dzx_rows)
              fprintf(fid,[indent,TD1,' = ',TD1,' + ',TD2,'(',Dind1,');\n']);
            else
              fprintf(fid,[indent,TD1,' = ',TD2,'(',Dind1,');\n']);
            end
          end
        else
          Dind2 = cadaindprint(dzyInds(:));
          fprintf(fid,[indent,TD1,'(',Dind2,') = ',...
            TD1,'(',Dind2,') + ',TD2,'(',Dind1,');\n']);
        end
      end
      fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
    end
  else
    nzlocs = [];
  end
elseif ~isempty(dxind)
  % ---------------------X has Derivatives----------------------------- %
  % we can compute DZ' = Y'*DX' where DX' is reshaped.
  
  % Get DZ from X.
  xTranMap = zeros(xMrow,xNcol); xTranMap(:) = 1:xMrow*xNcol;
  xTranMap = xTranMap.';
  dx       = sparse(dxind(:,1),dxind(:,2),1:nzx,xMrow*xNcol,nv);
  dxTran   = dx(xTranMap(:),:); % DX.'
  dxTranR  = reshape(dxTran,xNcol,xMrow*nv);
  switch caller
    case 'mtimes'
      dzxTranR = ytemp.'*dxTranR;
      dzxTran  = reshape(dzxTranR,FMrow*FNcol,nv);
    case 'mrdivide'
      dzxTranR = ytemp.'\dxTranR;
      dzxTran  = reshape(dzxTranR,FMrow*FNcol,nv);
    case 'mldivide'
      error('X should not have derivatives if called from mldivide')
  end
  
  zTranMap = zeros(FNcol,FMrow); zTranMap(:) = 1:FMrow*FNcol;
  zTranMap = zTranMap.';
  dzx      = dzxTran(zTranMap(:),:);
  [dzx_rows, dzx_cols] = find(dzx);
  if size(dzx_rows,2) > 1;dzx_rows = dzx_rows.';dzx_cols = dzx_cols.';end
  
  if ~isempty(dzx_rows)
    nzlocs = [dzx_rows,dzx_cols];
    if DPFLAG
      % ------------------Derivative Printing-------------------------- %
      % We need to project the derivatives of DX into a matrix of the
      % form of reshape(dxTran,xNcol,xMrow*nv) - but remove zeros columns
      TD1 = ['cada',NDstr,'td1'];
      
      % Determine how many columns of DX' reshaped actually have non-zero
      % derivatives
      dxTranRsum = sum(dxTranR,1);
      dxTranRsumlocs = find(dxTranRsum);
      coldim = length(dxTranRsumlocs);
      if coldim < xMrow*nv
        [dxTRrows, dxTRcols, dxTRlocs] = find(dxTranR(:,dxTranRsumlocs));
      else
        [dxTRrows, dxTRcols, dxTRlocs] = find(dxTranR); 
      end
      [~,dxTRlocs] = sortrows(dxTRlocs(:));
      dxTRrows = dxTRrows(dxTRlocs);
      dxTRcols = dxTRcols(dxTRlocs);
      
      % Project dx nonzeros into matrix size xN x coldim
      if xNcol*coldim > 250 && nzx < (3/4)*xNcol*coldim
        SPXflag = 1;
        % Do the projection sparsely - dxTran is larger than 250 and
        % has at least 25% zeros.
         TDind1 = cadaindprint(dxTRrows);
         TDind2 = cadaindprint(dxTRcols);
        % Print out the Projection
        fprintf(fid,[indent,TD1,' = sparse(',TDind1,',',TDind2,',',...
          x.deriv(Vcount).name,',%1.0d,%1.0d);\n'],xNcol,coldim);
      else
        SPXflag = 0;
        % Just project into zeros using a linear index - linear index
        % of dxTran and dxTran-reshaped are the same
        TDind1 = cadaindprint(sub2ind([xNcol,coldim],dxTRrows,dxTRcols));
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,%1.0d);\n'],...
          xNcol,coldim);
        fprintf(fid,[indent,TD1,'(',TDind1,') = ',...
          x.deriv(Vcount).name,';\n']);
      end
      
      % Can now print out Y.'*DX.' - where DX.' is in file as TD2.
      if strcmp(caller,'mrdivide')
        fprintf(fid,[indent,TD1,' = ',y.func.name,'.''\\',TD1,';\n']);
      else
        fprintf(fid,[indent,TD1,' = ',y.func.name,'.''*',TD1,';\n']);
      end
      
      % Need to get the mapping from dzxTran to dzx.
      if coldim < xMrow*nv
        [dzxTRrows,dzxTRcols] = find(dzxTranR(:,dxTranRsumlocs));
        dzxTranR = sparse(dzxTRrows,dxTranRsumlocs(dzxTRcols),...
          sub2ind([FNcol,coldim],dzxTRrows,dzxTRcols),FNcol,FMrow*nv);
        dzxTran = reshape(dzxTranR,FMrow*FNcol,nv);
      else
        [dzxTrows,dzxTcols] = find(dzxTran);
        dzxTran  = sparse(dzxTrows,dzxTcols,...
          sub2ind([FMrow*FNcol,nv],dzxTrows,dzxTcols),FMrow*FNcol,nv);
      end
      dzx      = dzxTran(zTranMap(:),:);
      dzxTinds = nonzeros(dzx);
      Dind1    = cadaindprint(dzxTinds(:));
      % dzxTinds are now the Linear Reference index off of dzxTran (in
      % the file) to the nonzeros in dzx.
      fprintf(fid,[indent,TD1,' = ',TD1,'(:);\n']);
      if SPXflag
        fprintf(fid,[indent,derivstr,' = full(',TD1,'(',Dind1,'));\n']);
      else
        fprintf(fid,[indent,derivstr,' = ',TD1,'(',Dind1,');\n']);
      end
    end
  else
    nzlocs = [];
  end
elseif ~isempty(dyind)
  % ---------------------Y has Derivatives----------------------------- %
  % We say DZ = X*DY - where DY is reshaped.
  % If mldivide: DZ = X\DY
  % Get DZ from Y first.
  dy  = sparse(dyind(:,1),dyind(:,2),1:nzy,yMrow*yNcol,nv);
  dyR = reshape(dy,yMrow,yNcol*nv);
  switch caller
    case 'mtimes'
      dzyR = xtemp*dyR;
      dzy  = reshape(dzyR,FMrow*FNcol,nv);
    case 'mldivide'
      dzyR = xtemp\dyR;
      dzy  = reshape(dzyR,FMrow*FNcol,nv);
    case 'mrdivide'
      error('Y should not have derivatives if called from mrdivide')
  end
  
  [dzy_rows, dzy_cols] = find(dzy);
  if size(dzy_rows,2) > 1;dzy_rows = dzy_rows.';dzy_cols = dzy_cols.';end
  if ~isempty(dzy_rows)
    nzlocs = [dzy_rows,dzy_cols];
    if DPFLAG
      % ------------------Derivative Printing-------------------------- %
      % Print out the Calculations for DZY
      TD1 = ['cada',NDstr,'td1'];
      % Need to project the vector of DY(in the file) into the matrix
      % of the form reshape(DY,yMrow,yNcol*nv) - except remove any known
      % zero columns
      
      % Determine how many columns of DY reshaped actually have non-zero
      % derivatives
      dyRsum = sum(dyR,1);
      dyRsumlocs = find(dyRsum);
      coldim = length(dyRsumlocs);

      if coldim < yNcol*nv
        % dyR has some columns which are all zero - do not need to sum
        % these columns
        dycaseflag = 'columns';
        [dyRrows, dyRcols] = find(dyR(:,dyRsumlocs));
      else
        % X is considered full, dyR has all non-zero columns, just project
        % and compute dzyR = X*dyR.
        dycaseflag = 'full';
        [dyRrows, dyRcols] = find(dyR);
      end
      
      % Project dy nonzeros into matrix size yM by coldim
      if yMrow*coldim > 250 && nzy < (3/4)*yMrow*coldim
        SPYflag = 1;
        % Do the projection sparsely - dy is larger than 250 and has at
        % least %25 zeros.
        TDind1 = cadaindprint(dyRrows);
        TDind2 = cadaindprint(dyRcols);
        % Print out the Projection
        fprintf(fid,[indent,TD1,' = sparse(',TDind1,',',TDind2,',',...
          y.deriv(Vcount).name,',%1.0d,%1.0d);\n'],yMrow,coldim);
      else
        SPYflag = 0;
        % Project into zeros using a linear index.
        TDind1 = cadaindprint(sub2ind([yMrow,coldim],dyRrows,dyRcols));
        % Print out the Projection
        fprintf(fid,[indent,TD1,' = zeros(%1.0d,%1.0d);\n'],...
          yMrow,coldim);
        fprintf(fid,[indent,TD1,'(',TDind1,') = ',...
          y.deriv(Vcount).name,';\n']);
      end
      
      % Can now print out X*DY - where DY is TD2
      if strcmp(caller,'mldivide')
        fprintf(fid,[indent,TD1,' = ',x.func.name,'\\',TD1,';\n']);
      else
        fprintf(fid,[indent,TD1,' = ',x.func.name,'*',TD1,';\n']);
      end
      % TD2 in the file now corresponds with DZY (except reshaped), so
      % we need to know what are the nonzero indices of DZY which we
      % need to reference off of TD2 (even though it is reshaped,
      % linear indexing will still be the same)
      if coldim < yNcol*nv
        dyInds = find(dzyR(:,dyRsumlocs));
      else
        dyInds = sub2ind([FMrow*FNcol,nv],dzy_rows,dzy_cols);
      end
      Dind1  = cadaindprint(dyInds(:));
      % dyInds are now the linear reference off of TD2 which give the
      % nonzeros in DZY
      fprintf(fid,[indent,TD1,' = ',TD1,'(:);\n']);
      if SPYflag
        fprintf(fid,[indent,derivstr,' = full(',TD1,'(',Dind1,'));\n']);
      else
        fprintf(fid,[indent,derivstr,' = ',TD1,'(',Dind1,');\n']);
      end
    end
  else
    nzlocs = [];
  end
end