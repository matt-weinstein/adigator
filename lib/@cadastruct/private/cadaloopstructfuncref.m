function cadaloopstructfuncref(inputstruct)
% This function prints out the iteration dependent structure/cell array
% referencing and assignment stuff - called from overloaded CADASTRUCT
% subsasgn/subsref
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

refflag = inputstruct.refflag;
% refflag = 1 => subsref => y = x(i)
% refflag = 0 => subsasgn => x(i) = y
funcstr    = inputstruct.yStr;
funrefstr  = inputstruct.xiStr;

fid        = inputstruct.fid;
indent     = inputstruct.indent;
sizecheck  = inputstruct.sizecheck;

if sizecheck
  IndName    = inputstruct.IndName;
  
  CountName  = inputstruct.CountName;
  vvec       = inputstruct.vvec;
  vMrow      = inputstruct.vsize(1);
  vNcol      = inputstruct.vsize(2);
  NDstr      = inputstruct.NDstr;
  vecDim     = inputstruct.vecDim;
  emptyfieldcheck = inputstruct.emptyfieldcheck;
  TF1 = ['cada',NDstr,'tempf1'];
  TF2 = ['cada',NDstr,'tempf2'];
  if ~strcmp(TF1,funrefstr)
    fprintf(fid,[indent,TF1,' = ',funrefstr,';\n']);
  end
  if refflag
    % Initialize to zeros if reference.
    if vvec == 1
      fprintf(fid,[indent,TF2,' = zeros(',vecDim,',%1.0f);\n'],vNcol);
    elseif vvec == 2
      fprintf(fid,[indent,TF2,' = zeros(%1.0f,',vecDim,');\n'],vMrow);
    else
      fprintf(fid,[indent,TF2,' = zeros(%1.0f,%1.0f);\n'],vMrow,vNcol);
    end
  end
  if emptyfieldcheck
    fprintf(fid,[indent,'cadaconditional1 = ~isempty(',TF1,');\n']);
    fprintf(fid,[indent,'if cadaconditional1\n']);
    indent = [indent,'    '];
  end
  if ~isempty(IndName)
    % Check to make sure that this reference actually happens
    IndDep = inputstruct.IndDep;
    if IndDep
      fprintf(fid,[indent,'cadaconditional1 = ~',IndName,'(%1.0f,',CountName,');\n'],inputstruct.Scount);
    else
      fprintf(fid,[indent,'cadaconditional1 = ~',IndName,'(%1.0f);\n'],inputstruct.Scount);
    end
    fprintf(fid,[indent,'if cadaconditional1\n']);
    indent = [indent,'    ']; 
  end
  if refflag
    if vvec == 1
      fprintf(fid,[indent,TF2,'(:,1:size(',TF1,',2)) = ',TF1,';\n']);
    elseif vvec == 2
      fprintf(fid,[indent,TF2,'(1:size(',TF1,',1),:) = ',TF1,';\n']);
    else
      fprintf(fid,[indent,TF2,'(1:size(',TF1,',1),1:size(',TF1,',2)) = ',TF1,';\n']);
    end
  else
    if vvec == 1
      fprintf(fid,[indent,TF2,' = ',funcstr,'(:,1:size(',TF1,',2));\n']);
    elseif vvec == 2
      fprintf(fid,[indent,TF2,' = ',funcstr,'(1:size(',TF1,',1),:);\n']);
    else
      fprintf(fid,[indent,TF2,' = ',funcstr,'(1:size(',TF1,',1),1:size(',TF1,',2));\n']);
    end
  end
  if ~isempty(IndName)
    indent(1:4) = [];
    fprintf(fid,[indent,'end\n']);
  end
  if emptyfieldcheck
    indent(1:4) = [];
    fprintf(fid,[indent,'end\n']);
  end
  if refflag
    fprintf(fid,[indent,funcstr,' = ',TF2,';\n']);
  else
    fprintf(fid,[indent,funrefstr,' = ',TF2,';\n']);
  end
elseif refflag
  fprintf(fid,[indent,funcstr,' = ',funrefstr,';\n']);
else
  fprintf(fid,[indent,funrefstr,' = ',funcstr,';\n']);
end