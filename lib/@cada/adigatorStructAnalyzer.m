function x = adigatorStructAnalyzer(x,xStr,subsflag)
% CADA overloaded version of adigatorStructAnalyzer
%
% Other versions are called to recursively parse structure/cell objects and
% call adigatorVarAnalyzer on each of the CADA objects. The fact that this
% has been called implies that the CADA object lives within a
% cell/structure, thus there are a few special cases that are checked for
% prior to callling adigatorVarAnalyzer.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
NUMvod  =   ADIGATOR.NVAROFDIFF;
xID = x.id;
PFLAG = ADIGATOR.RUNFLAG == 2 && ADIGATOR.PRINT.FLAG && ...
  ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,2) ~= ADIGATOR.VARINFO.NAMELOCS(xID,2);
if PFLAG
  fid = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
end
% Theres a special case where this needs to print stuff when calling a
% previously created derivative file
if size(ADIGATOR.VARINFO.NAMELOCS,1) >= xID 
  if isinf(ADIGATOR.VARINFO.NAMELOCS(xID,3))
    ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = ADIGATOR.VARINFO.NAMELOCS(xID,3);
  end
  if ADIGATOR.VARINFO.NAMELOCS(xID,2) < 0
    ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,2) = ADIGATOR.VARINFO.NAMELOCS(xID,2);
  end
end

ADIGATOR.VARINFO.LASTOCC([xID,ADIGATOR.VARINFO.COUNT],1) = ADIGATOR.VARINFO.COUNT;
if ADIGATOR.RUNFLAG > 0
  if isinf(ADIGATOR.VARINFO.NAMELOCS(xID,3))
    ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = ADIGATOR.VARINFO.NAMELOCS(xID,3);
  end
  if ADIGATOR.VARINFO.NAMELOCS(xID,2) < 0
    ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,2) = ADIGATOR.VARINFO.NAMELOCS(xID,2);
  end
  [funcstr,DPFLAG] = cadafuncname();
  for Vcount = 1:NUMvod
    if ~isempty(x.deriv(Vcount).nzlocs)
      derivstr = cadadername(funcstr,Vcount);
      if PFLAG && DPFLAG
        fprintf(fid,[indent,derivstr,' = ',x.deriv(Vcount).name,';\n']);
      end
      x.deriv(Vcount).name = derivstr;
    end
  end
  if PFLAG
    fprintf(fid,[indent,funcstr,' = ',x.func.name,';\n']);
  end
  x.func.name = funcstr;
end
x.id = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
x = adigatorVarAnalyzer('',x,xStr,subsflag);
if isinf(ADIGATOR.VARINFO.NAMELOCS(xID,3))
  ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT-1,3) = ADIGATOR.VARINFO.NAMELOCS(xID,3);
end
if ADIGATOR.VARINFO.NAMELOCS(xID,2) < 0
  ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT-1,2) = ADIGATOR.VARINFO.NAMELOCS(xID,2);
end
end