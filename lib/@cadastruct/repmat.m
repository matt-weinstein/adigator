function y = repmat(x,varargin)
% CADASTRUCT overloaded version of REPMAT
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;

% ------------------Parse Sizing Inputs---------------------------------- %
if nargin == 1
  error('Requires at least 2 inputs.')
elseif nargin == 2
  if isa(varargin{1},'cada')
    if varargin{1}.func.size(1)*varargin{1}.func.size(2) == 1
      if ~isempty(varargin{1}.func.value)
        repMrow = varargin{1}.func.value;
        repNcol = repMrow;
      elseif ~ADIGATOR.EMPTYFLAG
        error('Cannot REPMAT a purely symbolic dimension.')
      else
        repMrow = [];
        repNcol = [];
      end
    elseif varargin{1}.func.size(1)*varargin{1}.func.size(2) == 2
      if ~isempty(varargin{1}.func.value)
        repMrow = varargin{1}.func.value(1);
        repNcol = varargin{1}.func.value(2);
      elseif ~ADIGATOR.EMTPYFLAG
        error('Cannot REPMAT a purely symbolic dimension.')
      else
        repMrow = [];
        repNcol = [];
      end
    else
      error('Can only work with 2 dimensions')
    end
    ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(varargin{1})
    if length(varargin{1}) == 1
      repMrow = varargin{1};
      repNcol = repMrow;
    elseif length(varargin{1}) == 2
      repMrow = varargin{1}(1);
      repNcol = varargin{1}(2);
    else
      error('Can only work with 2 dimensions')
    end
  else
    error('Invalid input');
  end
elseif nargin == 3
  if isa(varargin{1},'cada')
    if varargin{1}.func.size(1)*varargin{1}.func.size(2) == 1
      if ~isempty(varargin{1}.func.value)
        repMrow = varargin{1}.func.value;
      elseif ~ADIGATOR.EMPTYFLAG
        error('Cannot REPMAT a purely symbolic dimension.')
      else
        repMrow = [];
      end
    else
      error('Can only work with 2 dimensions')
    end
    ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(varargin{1})
    if length(varargin{1}) == 1
      repMrow = varargin{1};
    else
      error('Can only work with 2 dimensions')
    end
  else
    error('Invalid input');
  end
  if isa(varargin{2},'cada')
    if varargin{2}.func.size(1)*varargin{2}.func.size(2) == 1
      if ~isempty(varargin{2}.func.value)
        repNcol = varargin{2}.func.value;
      elseif ~ADIGATOR.EMPTYFLAG
        error('Cannot REPMAT a purely symbolic dimension.')
      else
        repNcol = [];
      end
    else
      error('Can only work with 2 dimensions')
    end
    ADIGATOR.VARINFO.LASTOCC(varargin{2}.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isnumeric(varargin{2})
    if length(varargin{2}) == 1
      repNcol = varargin{2};
    else
      error('Can only work with 2 dimensions')
    end
  else
    error('Invalid input');
  end
else
  error('Too many input arguments.')
end

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
if ~ADIGATOR.EMPTYFLAG || (~isempty(repMrow) && ~isempty(repNcol))
  % This likes to call CADA operations for some strange reason 
  precount = ADIGATOR.VARINFO.COUNT;
  preflag  = ADIGATOR.PRINT.FLAG;
  ADIGATOR.PRINT.FLAG = 0;
  y.val = repmat(x.val,repMrow,repNcol);
  ADIGATOR.PRINT.FLAG = preflag;
  ADIGATOR.VARINFO.COUNT = precount;
end


if ~x.arrayflag
  fnames = fieldnames(x.val);
  for Fcount = 1:length(fnames)
    F = fnames{Fcount};
    v = x.val.(F);
    if isa(v,'cada')
      ADIGATOR.VARINFO.LASTOCC(v.id,1) = ADIGATOR.VARINFO.COUNT;
      if ADIGATOR.RUNFLAG > 0 && isinf(ADIGATOR.VARINFO.NAMELOCS(v.id,3))
        ADIGATOR.VARINFO.NAMELOCS(y.id,3) = ADIGATOR.VARINFO.NAMELOCS(v.id,3);
      end
    end
  end
else
  if ADIGATOR.RUNFLAG > 0 && isinf(ADIGATOR.VARINFO.NAMELOCS(x.id,3))
    ADIGATOR.VARINFO.NAMELOCS(y.id,3) = ADIGATOR.VARINFO.NAMELOCS(x.id,3);
  end
end
if PFLAG && ~ADIGATOR.EMPTYFLAG
  fprintf(fid,[indent,yname,' = repmat(',x.name,',%1.0d,%1.0d);\n'],repMrow,repNcol);
end

y.name = yname;
y.arrayflag = 1;

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
end