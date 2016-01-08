function adigatorPrintStructAsgn(x,xStrA,xStrR,idA,idR)
% this does a recursive call to build proper strings and print out
% structure assignments
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if isa(x,'cadastruct')
  x = x.val;
end

if isa(x,'cada')
  fid = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  namelocA = ADIGATOR.VARINFO.NAMELOCS(idA,1);
  Aflag = 0;
  if namelocA == 0
    Aflag = 1;
    namelocA = length(ADIGATOR.VARINFO.NAMES)+1;
    ADIGATOR.VARINFO.NAMES{namelocA} = [];
  end
  namelocR = ADIGATOR.VARINFO.NAMELOCS(idR,1);
  Rflag = 0;
  if namelocR == 0
    Rflag = 1;
    namelocR = length(ADIGATOR.VARINFO.NAMES)+1;
    ADIGATOR.VARINFO.NAMES{namelocR} = [];
  end
  oldnameA = ADIGATOR.VARINFO.NAMES{namelocA};
  oldnameR = ADIGATOR.VARINFO.NAMES{namelocR};
  ADIGATOR.VARINFO.NAMES{namelocA} = xStrA;
  ADIGATOR.VARINFO.NAMES{namelocR} = xStrR;
  [funcstrA,DPflag] = cadafuncname(idA);
  funcstrR          = cadafuncname(idR);
  xderiv = x.deriv;
  if DPflag
    for Vcount = 1:length(xderiv)
      if ~isempty(xderiv(Vcount).nzlocs)
        derivstrA = cadadername(funcstrA,Vcount,idA);
        derivstrR = cadadername(funcstrR,Vcount,idR);
        fprintf(fid,[indent,derivstrA,' = ',derivstrR,';\n']);
      end
    end
  end
  fprintf(fid,[indent,funcstrA,' = ',funcstrR,';\n']);
  ADIGATOR.VARINFO.NAMES{namelocA} = oldnameA;
  ADIGATOR.VARINFO.NAMES{namelocR} = oldnameR;
  if Rflag
    ADIGATOR.VARINFO.NAMES(namelocR) = [];
  end
  if Aflag
    ADIGATOR.VARINFO.NAMES(namelocA) = [];
  end
elseif isstruct(x)
  fnames = fieldnames(x);
  for Fcount = 1:length(fnames)
    F = fnames{Fcount};
    if numel(x) > 1
      for I = 1:size(x,1)
        for J = 1:size(x,2)
          Istr = sprintf('%1.0d,%1.0d',I,J);
          adigatorPrintStructAsgn(x(I,J).(F),[xStrA,'(',Istr,').',F],[xStrR,'(',Istr,').',F],idA,idR);
        end
      end
    else
      adigatorPrintStructAsgn(x.(F),[xStrA,'.',F],[xStrR,'.',F],idA,idR);
    end
  end
elseif iscell(x)
  for I = 1:size(x,1)
    for J = 1:size(x,2)
      Istr = sprintf('%1.0d,%1.0d',I,J);
      adigatorPrintStructAsgn(x{I,J},[xStrA,'{',Istr,'}'],[xStrR,'{',Istr,'}'],idA,idR);
    end
  end
else
  fid = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  fprintf(fid,[indent,xStrA,' = ',xStrR,';\n']);
end

end