function adigatorPrintOutputIndices(x)
% Parses through object to call CADA version of adigatorPrintOutputIndices
%
% See also: cada/adigatorPrintOutputIndices
ParsePrint(x.id,x.val,x.name);
end

function ParsePrint(sid,x,xname)
global ADIGATOR
if isa(x,'cada')
  % Need to give this the proper names - going to do this by creating
  % ''dummy'' entries in NAMELOCS
  xid = ADIGATOR.VARINFO.COUNT;
  adigatorAssignImpVarNames(xid,xname,0);
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  x.id = xid;
  ADIGATOR.VARINFO.NAMELOCS(xid,2) = ADIGATOR.VARINFO.NAMELOCS(sid,2);
  adigatorPrintOutputIndices(x);
elseif isstruct(x) && ~isempty(x)
  fnames = fieldnames(x);
  if numel(x) > 1
    for I = 1:numel(x)
      for J = 1:length(fnames)
        F = fnames{J};
        ParsePrint(sid,x(I).(F),sprintf([xname,'(%1.0f).',F],I));
      end
    end
  else
    for J = 1:length(fnames)
      F = fnames{J};
      ParsePrint(sid,x.(F),[xname,'.',F]);
    end
  end
elseif iscell(x) && ~isempty(x)
  for I = 1:numel(x)
    ParsePrint(sid,x{I},sprintf([xname,'{%1.0f}'],I));
  end
end
end