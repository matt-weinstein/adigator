function y = ctranspose(x)
% CADASTRUCT overloaded version of CTRANSPOSE
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;

y = x;
y.id = ADIGATOR.VARINFO.COUNT;
if ADIGATOR.RUNFLAG == 2
  nameloc = ADIGATOR.VARINFO.NAMELOCS(y.id,1);
  if nameloc > 0
    yname = ADIGATOR.VARINFO.NAMES{nameloc};
  else
    yname = sprintf(['cada',NDstr,'s%1.0f'],ADIGATOR.VARINFO.NAMELOCS(yid,2));
  end
else
  yname = 'cadadummystruct';
end
y.name = yname;
y.val = ctranspose(x.val);

if PFLAG && ~ADIGATOR.EMPTYFLAG
  fprintf(fid,[indent,yname,' = ctranspose(',x.name,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
end