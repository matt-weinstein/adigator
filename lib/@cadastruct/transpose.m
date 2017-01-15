function y = transpose(x)
% CADASTRUCT overloaded version of TRANSPOSE
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
    yname = sprintf('cada%1.0ds%1.0f',ADIGATOR.NVAROFDIFF,ADIGATOR.VARINFO.NAMELOCS(y.id,2));
  end
else
  yname = 'cadadummystruct';
end
y.name = yname;
y.val = transpose(x.val);

if PFLAG && ~ADIGATOR.EMPTYFLAG
  fprintf(fid,[indent,yname,' = transpose(',x.name,');\n']);
end
if ADIGATOR.RUNFLAG > 0 && isinf(ADIGATOR.VARINFO.NAMELOCS(x.id,3))
  ADIGATOR.VARINFO.NAMELOCS(y.id,3) = ADIGATOR.VARINFO.NAMELOCS(x.id,3);
end
ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
end