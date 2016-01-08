function cadaloopstructderivref(inputstruct)
% This function prints out the iteration dependent structure/cell array
% referencing and assignment stuff - called from overloaded CADASTRUCT
% subsasgn/subsref
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

refflag = inputstruct.refflag;
% refflag = 1 => subsref => y = x(i)
% refflag = 0 => subsasgn => x(i) = y
if refflag
  derivstr    = inputstruct.dyStr;
  derivrefstr = inputstruct.dxiStr;
else
  derivstr    = inputstruct.dxiStr;
  derivrefstr = inputstruct.dyStr;
end

fid        = inputstruct.fid;
indent     = inputstruct.indent;
IndName    = inputstruct.IndName;
IndFlags   = inputstruct.IndFlags;
CountName  = inputstruct.CountName;
vvec       = inputstruct.vvec;
nzd        = inputstruct.nzd;
NDstr      = inputstruct.NDstr;
vecDim     = inputstruct.vecDim;
emptycheck = inputstruct.emptycheck;
emptyfieldcheck = inputstruct.emptyfieldcheck; 
% this will only be true if being called from s(i).y(j) = b, where b is cada

if isempty(IndName)
  fprintf(fid,[indent,derivstr,' = ',derivrefstr,';\n']);
else
  if vvec
    % Vectorized
    if IndFlags(1)
      % Indices change on this loop
      asgnind = ['(:,logical(',IndName,'(:,',CountName,')))'];
    else
      % Indices do not change on this loop
      asgnind = ['(:,logical(',IndName,'))'];
    end
  else
    % Non-Vectorized
    if IndFlags(1)
      % Indices change on this loop
      asgnind = ['(logical(',IndName,'(:,',CountName,')))'];
    else
      % Indices do not change on this loop
      asgnind = ['(logical(',IndName,'))'];
    end
  end
  
  if refflag
    % Initialize dy if reference - dont need initiliaze if is assignment
    TD1 = ['cada',NDstr,'td1'];
    if vvec
      fprintf(fid,[indent,TD1,' = zeros(',vecDim,',%1.0d);\n'],nzd);
    else
      fprintf(fid,[indent,TD1,' = zeros(%1.0d,1);\n'],nzd);
    end
  end
  if emptycheck
    TF1 = ['cada',NDstr,'tf1'];
    if vvec
      fprintf(fid,[indent,TF1,' = ',asgnind(4:end-1),';\n']);
      asgnind = ['(:,',TF1,')'];
    else
      fprintf(fid,[indent,TF1,' = ',asgnind(2:end-1),';\n']);
      asgnind = ['(',TF1,')'];
    end
    fprintf(fid,[indent,'cadaconditional1 = any(',TF1,');\n']);
    fprintf(fid,[indent,'if cadaconditional1\n']);
    indent = [indent,'    '];
  elseif emptyfieldcheck
    TD2 = ['cada',NDstr,'td2'];
    fprintf(fid,[indent,TD2,' = ',derivrefstr,';\n']);
    fprintf(fid,[indent,'cadaconditional1 = ~isempty(',TD2,');\n']);
    fprintf(fid,[indent,'if cadaconditional1\n']);
    indent = [indent,'    '];
    derivrefstr = TD2;
  end
  if refflag
    % reference is dy(logical) = dxi
    fprintf(fid,[indent,TD1,asgnind,' = ',derivrefstr,';\n']);
  else
    % assignment is dxi = dy(logical)
    fprintf(fid,[indent,derivstr,' = ',derivrefstr,asgnind,';\n']);
  end
  if emptycheck || emptyfieldcheck
    indent(1:4) = [];
    fprintf(fid,[indent,'end\n']);
  end
  if refflag
    fprintf(fid,[indent,derivstr,' = ',TD1,';\n']);
  end
end