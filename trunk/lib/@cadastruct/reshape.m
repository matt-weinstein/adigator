function y = reshape(x,varargin)
% CADASTRUCT overloaded version of function RESHAPE.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
fid    = ADIGATOR.PRINT.FID;
PFLAG  = ADIGATOR.PRINT.FLAG;
indent = ADIGATOR.PRINT.INDENT;
% ----------------------------Determine Inputs----------------------------%
varargSize = length(varargin);
if varargSize == 0 ;
  error('Requires at least 2 inputs.');
elseif varargSize == 1
  if isa(varargin{1},'cada')
    if ~isempty(varargin{1}.func.value)
      ADIGATOR.VARINFO.LASTOCC(varargin{1}.id,1) = ADIGATOR.VARINFO.COUNT;
      if length(varargin{1}.func.value) > 2
        error('cannot RESHAPE more than 2 dimensions')
      elseif length(varargin{1}.func.value) == 1
        error('Size vector must have at least two elements.');
      end
      FMrow = varargin{1}.func.value(1);
      FNcol = varargin{1}.func.value(2);
      RepStr = varargin{1}.func.name;
    elseif ~ADIGATOR.EMPTYFLAG
      error('cannot RESHAPE a strictly symbolic dimension')
    end
  elseif length(varargin{1}) == 2
    FMrow  = varargin{1}(1);
    FNcol  = varargin{1}(2);
    RepStr = sprintf('%1.0f,%1.0f',FMrow,FNcol);
  elseif ~ADIGATOR.EMPTYFLAG
    error('cannot RESHAPE more than 2 dimensions')
  else
    FMrow = [];
    FNcol = [];
  end
elseif varargSize == 2
  FMrow = varargin{1};
  if isa(FMrow,'cada')
    if ~isempty(FMrow.func.value)
      ADIGATOR.VARINFO.LASTOCC(FMrow.id,1) = ADIGATOR.VARINFO.COUNT;
      RowStr = FMrow.func.name;
      FMrow  = FMrow.func.value;
    elseif ~ADIGATOR.EMPTYFLAG
      error('cannot RESHAPE a strictly symbolic dimension')
    else
      FMrow = [];
    end
  else
    RowStr = sprintf('%1.0f',FMrow);
  end
  FNcol = varargin{2};
  if isa(FNcol,'cada')
    if ~isempty(FNcol.func.value)
      ADIGATOR.VARINFO.LASTOCC(FNcol.id,1) = ADIGATOR.VARINFO.COUNT;
      ColStr = FNcol.func.name;
      FNcol  = FNcol.func.value;
    elseif ~ADIGATOR.EMPTYFLAG
      error('cannot RESHAPE a strictly symbolic dimension')
    else
      FNcol = [];
    end
  else
    ColStr = sprintf('%1.0f',FNcol);
  end
  RepStr = ['',RowStr,',',ColStr,''];
else
  error('Too many input arguments.');
end
if FMrow*FNcol ~= numel(x.val)
  error('To RESHAPE the number of elements must not change.')
end

% --------------------------Parse X-------------------------------------- %
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
if ~ADIGATOR.EMPTYFLAG || (~isempty(FMrow) && ~isempty(FNcol))
  y.val = reshape(x.val,FMrow,FNcol);
end

if ADIGATOR.RUNFLAG > 0 && isinf(ADIGATOR.VARINFO.NAMELOCS(x.id,3))
  ADIGATOR.VARINFO.NAMELOCS(y.id,3) = ADIGATOR.VARINFO.NAMELOCS(x.id,3);
end
if PFLAG && ~ADIGATOR.EMPTYFLAG
  fprintf(fid,[indent,yname,' = reshape(',x.name,',',RepStr,');\n']);
end

ADIGATOR.VARINFO.LASTOCC([y.id x.id],1) = ADIGATOR.VARINFO.COUNT;
ADIGATOR.VARINFO.COUNT                  = ADIGATOR.VARINFO.COUNT+1;
end