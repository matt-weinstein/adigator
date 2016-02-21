function y = horzcat(varargin)
% CADASTRUCT overloaded version of function HORZCAT.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
numder = ADIGATOR.DERNUMBER;
X    = cell(1,nargin);
Xstr = cell(1,nargin);
% See which are cadastructs, convert those which are not.
TVcount = 0;
for I = 1:nargin
  xi = varargin{I};
  if isa(xi,'cadastruct')
    X{I} = xi.val;
    Xstr{I} = [xi.name,','];
    ADIGATOR.VARINFO.LASTOCC(xi.id,1) = ADIGATOR.VARINFO.COUNT;
  else
    TVcount = TVcount+1;
    xiStr   = sprintf('cada%1.0dtemps%1.d',numder,TVcount);
    xi      = parseinput(xi,xiStr);
    X{I}    = xi;
    Xstr{I} = [xiStr,','];
  end
end
Xstr{end}(end) = [];

y.id = ADIGATOR.VARINFO.COUNT;
%y.id = ADIGATOR.VARINFO.COUNT;
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
if ADIGATOR.EMPTYFLAG
  try
    y.val = horzcat(X{:});
  catch
    y.val = X{1};
  end
else
  y.val = horzcat(X{:});
end
y.arrayflag = 1;

if PFLAG && ~ADIGATOR.EMPTYFLAG
  fprintf(fid,[indent,yname,' = horzcat(',Xstr{:},');\n']);
end
ADIGATOR.VARINFO.LASTOCC(y.id,1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT           = ADIGATOR.VARINFO.COUNT+1;
y = cadastruct(y);
end

function yi = parseinput(xi,xiStr)
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
switch class(xi)
  case 'cada'
    ADIGATOR.VARINFO.LASTOCC(xi.id,1) = ADIGATOR.VARINFO.COUNT;
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      yid = ADIGATOR.VARINFO.COUNT;
      ynameloc = ADIGATOR.VARINFO.NAMELOCS(yid,1);
      if ynameloc > 0
        nameloc = ynameloc;
        oldname = ADIGATOR.VARINFO.NAMES{nameloc};
      else
        oldnameinfo = ADIGATOR.VARINFO.NAMELOCS(yid,:);
        nameloc     = length(ADIGATOR.VARINFO.NAMES);
        oldname     = ADIGATOR.VARINFO.NAMES{nameloc};
        newnameinfo = [nameloc 0 1];
        ADIGATOR.VARINFO.NAMELOCS(yid,:) = newnameinfo;
      end
      ADIGATOR.VARINFO.NAMES{nameloc} = xiStr;
      
      [funcstr,DPflag] = cadafuncname(yid);
      if DPflag
        xideriv = xi.deriv;
        for Vcount = 1:length(xideriv)
          if ~isempty(xideriv(Vcount).nzlocs)
            derivstr = cadadername(funcstr,Vcount,yid);
            fprintf(fid,[indent,derivstr,' = ',xideriv(Vcount).name,';\n']);
          end
        end
      end
      
      fprintf(fid,[indent,funcstr,' = ',xi.func.name,';\n']);

      ADIGATOR.VARINFO.NAMES{nameloc} = oldname;
      if ynameloc == 0;
        ADIGATOR.VARINFO.NAMELOCS(yid,:) = oldnameinfo;
      end
    end
    yi = xi;
  case 'cadastruct'
    ADIGATOR.VARINFO.LASTOCC(xi.id,1) = ADIGATOR.VARINFO.COUNT;
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,xiStr,' = ',xi.name,';\n']);
    end
    yi = xi.val;
  case 'struct'
    fnames = fieldnames(xi);
    xiSize = size(xi);
    if length(xiSize) > 2
      error('ADiGator can only use 2D cells/structures')
    end
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      structinit = repmat({sprintf(',cell(%1.0f,%1.0f),',xiSize)},1,length(fnames));
      structinit{end}(end) = [];
      structinit = [fnames.';structinit];
      for F = 1:length(fnames)
        structinit{1,F} = ['''',structinit{1,F},''''];
      end
      fprintf(fid,[indent,xiStr,' = struct(',structinit{:},');\n']);
    end
    yi = xi;
    for Fcount = 1:length(fnames)
      F = fnames{Fcount};
      if numel(xi) == 1
        yi.(F) = parseinput(xi.(F),[xiStr,'.',F]);
      else
        for I = 1:numel(xi)
          yi(I).(F) = parseinput(xi(I).(F),sprintf([xiStr,'(%1.0d).',F],I));
        end
      end
    end
  case 'cell'
    xiSize = size(xi);
    if length(xiSize) > 2
      error('ADiGator can only use 2D cells/structures')
    end
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,xiStr,' = cell(%1.0f,%1.0f);\n'],xiSize);
    end
    yi = xi;
    for I = 1:xiSize(1)*xiSize(2)
      yi{I} = parseinput(xi{I},sprintf([xiStr,'{%1.0d}'],I));
    end
  case 'double'
    if PFLAG && ~ADIGATOR.EMPTYFLAG && ~isempty(xi)
      if numel(xi) < 10 && all(floor(xi(:)) == xi(:))
        xiname = mat2str(b);
      else
        xiname = cadamatprint(b);
      end
      fprintf(fid,[indent,xiStr,' = ',xiname,';\n']);
    end
    yi = xi;
  case 'char'
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,xiStr,' = ''',xi,''';\n']);
    end
    yi = xi;
end  
end
