function x = adigatorStructAnalyzer(x,xStr,subsflag)
% CADASTRUCT overloaded version of adigatorStructAnalyzer
%
% Called to recursively parse structure/cell objects and call
% adigatorVarAnalyzer on each of the CADA objects.

global ADIGATOR
if isa(x,'cadastruct') && ~x.arrayflag
  x = x.val;
  if ~isempty(x)
    x = orderfields(x);
    fnames = fieldnames(x);
    for I = 1:length(fnames)
      F = fnames{I};
      x.(F) = adigatorStructAnalyzer(x.(F),[xStr,'.',F],0);
    end
  end
else
  xID = x.id;
  if size(ADIGATOR.VARINFO.NAMELOCS,1) >= xID
    if isinf(ADIGATOR.VARINFO.NAMELOCS(xID,3))
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = ADIGATOR.VARINFO.NAMELOCS(xID,3);
    end
    if ADIGATOR.VARINFO.NAMELOCS(xID,2) < 0
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,2) = ADIGATOR.VARINFO.NAMELOCS(xID,2);
    end
  end
  ADIGATOR.VARINFO.LASTOCC([xID,ADIGATOR.VARINFO.COUNT],1) = ADIGATOR.VARINFO.COUNT;
  xid = ADIGATOR.VARINFO.COUNT;
  if ADIGATOR.RUNFLAG == 2 && ~ADIGATOR.EMPTYFLAG && ...
      any(ADIGATOR.VARINFO.NAMELOCS(xid,2) ~= ...
      ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.LASTOCC(:,1)==xid,2))
    adigatorPrintStructAsgn(x,xStr,x.name,xid,x.id);
  end
  x.name = xStr;
  x.id   = xid;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
  x = adigatorVarAnalyzer('',x,xStr,subsflag);
end
