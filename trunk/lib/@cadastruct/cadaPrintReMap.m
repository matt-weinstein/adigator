function y = cadaPrintReMap(x,y,varID)
% CADASTRUCT version of cadaPrintReMap
%
% Recursively looks through object to find CADA variables and call CADA
% version of cadaPrintReMap.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR


y.id = varID;
y.name = ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(varID,1)};
remapdriver(x,y,y.name,varID);
end

function remapdriver(x,y,name,varID)
global ADIGATOR
xarrayflag = 0;
yarrayflag = 0;
if isa(x,'cadastruct')
  xarrayflag = x.arrayflag;
  x = x.val;
end
if isa(y,'cadastruct')
  yarrayflag = y.arrayflag;
  y = y.val;
end

xclass = class(x); yclass = class(y);
if ~strcmp(xclass,yclass)
  if (isa(x,'cada') && prod(x.func.size)  == 0) || (isnumeric(x) && isempty(x))
    switch yclass
      case 'struct'
        x = struct();
      case 'cell'
        x = {};
      otherwise
       x = [];
    end
  elseif (isa(y,'cada') && prod(y.func.size)  == 0) || (isnumeric(y) && isempty(y))
    switch xclass
      case 'struct'
        y = struct();
      case 'cell'
        y = {};
      otherwise
       y = [];
    end
    yclass = xclass;
  else
    error([name,' takes on either an ',xclass,' or a ',yclass,...
      ' depending upon a loop iteration, conditional statement,',...
      ' function call, etc. This is not allowed']);
  end
elseif isa(x,'cada') && all(x.func.size == 0)
  x = [];
elseif isa(y,'cada') && all(y.func.size == 0)
  y = [];
end
switch yclass
  case 'cada'
    % change NAMES up
    nameloc = ADIGATOR.VARINFO.NAMELOCS(varID,1);
    oldname = ADIGATOR.VARINFO.NAMES{nameloc};
    ADIGATOR.VARINFO.NAMES{nameloc} = name;
    cadaPrintReMap(x,y,varID);
    ADIGATOR.VARINFO.NAMES{nameloc} = oldname;
  case 'struct'
    structureremap(x,y,name,varID,or(xarrayflag,yarrayflag));
  case 'cell'
    cellremap(x,y,name,varID);
end


end

function structureremap(x,y,name,varID,arrayflag)
global ADIGATOR
fid = ADIGATOR.PRINT.FID; indent = ADIGATOR.PRINT.INDENT;

xsize = size(x); xN = numel(x); xfields = fieldnames(x);
ysize = size(y); yN = numel(y); yfields = fieldnames(y);

[xyfields,xi,yi] = intersect(xfields,yfields);
xfields(xi) = [];
yfields(yi) = [];
% xyfields - common fields
% xfields  - fields only x has
% yfields  - fields only y has
% - due to union procedure, it should be the case that either xfields is
% empty or yfields is empty
% xfields not empty => y < x
% yfields not empty => x < y
if ~isempty(xfields) && ~isempty(yfields)
  error('Remapping error - please report to sourceforge forums')
end

csize = min(xsize,ysize);

if yN > 1 || xN > 1 || arrayflag
  % ------------------------ structure array ---------------------------- %
  
  % Take care of uncommon fields
  if ~isempty(yfields)
    % must add fields to x
    
    %ijstr = sprintf('(%1.0d,%1.0d)',ysize);
    for K = 1:length(yfields)
      F = yfields{K};
      %fprintf(fid,[indent,name,ijstr,'.',F,' = [];']);
      % Initialize any embedded stuff
      for I = ysize(1):-1:1
        for J = ysize(2):-1:1
          remapdriver([],y(I,J).(F),sprintf([name,'(%1.0d,%1.0d).',F],I,J),varID);
        end
      end
    end
    % Need to get anything in common fields too if larger size
    if any(ysize > xsize)
      for K = 1:length(xyfields)
        F = xyfields{K};
        % append rows
        for I = xsize(1)+1:ysize(1)
          for J = 1:xsize(2)
            remapdriver([],y(I,J).(F),sprintf([name,'(%1.0d,%1.0d).',F],I,J),varID);
          end
        end
        % append cols
        for I = 1:ysize(1)
          for J = xsize(2)+1:ysize(2)
            remapdriver([],y(I,J).(F),sprintf([name,'(%1.0d,%1.0d).',F],I,J),varID);
          end
        end
      end
    end
  elseif ~isempty(xfields)
    % must remove fields from x
    
    for I = 1:length(xfields)
      F = xfields{I};
      fprintf(fid,[indent,name,' = rmfield(',name,',''',F,''');\n']);
    end
  elseif any(ysize > xsize)
    % must add entries to x array
    
    %ijstr = sprintf('(%1.0d,%1.0d)',ysize);
    if ~isempty(xyfields)
      %F = xyfields{end};
      %fprintf(fid,[indent,name,ijstr,'.',F,' = [];\n']);
      % Initiliaze any embedded stuff
      for K = 1:length(xyfields)
         F = xyfields{K};
        % append rows
        for I = ysize(1):-1:xsize(1)+1
          for J = xsize(2):-1:1
            remapdriver([],y(I,J).(F),sprintf([name,'(%1.0d,%1.0d).',F],I,J),varID);
          end
        end
        % append cols
        for I = ysize(1):-1:1
          for J = ysize(2):-1:xsize(2)+1
            remapdriver([],y(I,J).(F),sprintf([name,'(%1.0d,%1.0d).',F],I,J),varID);
          end
        end
      end
    else
      fprintf(fid,[indent,name,ijstr,' = struct();\n']);
    end
  elseif any(xsize > ysize)
    % must remove entries from x array
    
    fprintf(fid,[indent,name,' = ',name,'(1:%1.0f,1:%1.0f);\n'],ysize);
  end
  
  % Remap common entries
  for K = 1:length(xyfields)
    F = xyfields{K};
    for I = csize(1):-1:1
      for J = csize(2):-1:1
        remapdriver(x(I,J).(F),y(I,J).(F),sprintf([name,'(%1.0d,%1.0d).',F],I,J),varID);
      end
    end
  end
elseif yN == 1
  % -------------------------- Scalar structure ------------------------- %
  if ~isempty(yfields)
    % must add fields to x
    for I = 1:length(yfields)
      F = yfields{I};
      % Don't print anything out here - let it get printed later
      remapdriver([],y.(F),[name,'.',F],varID);
    end
  elseif ~isempty(xfields)
    % must remove fields from x
    for I = 1:length(xfields)
      F = xfields{I};
      fprintf(fid,[indent,name,' = rmfield(',name,',''',F,''');\n']);
    end
  end
  % Remap common entries
  for K = 1:length(xyfields)
    F = xyfields{K};
    if xN == 1
      remapdriver(x.(F),y.(F),[name,'.',F],varID);
    else
      remapdriver([],y.(F),[name,'.',F],varID);
    end
  end
elseif yN == 0
  fprintf(fid,[indent,name,' = struct();\n']);
end
end

function cellremap(x,y,name,varID)
xsize = size(x); ysize = size(y);

if any(ysize > xsize)
  % Adding cell elements
  csize = xsize;
  % Initialize Cell Size
  %fprintf(fid,[indent,name,'{%1.0f,%1.0f} = [];\n'],ysize);

  % Append Rows
  for I = ysize(1):-1:xsize(1)+1
    for J = ysize(2):-1:1
      remapdriver([],y{I,J},sprintf([name,'{%1.0f,%1.0f}'],I,J),varID);
    end
  end
  % Append Columns
  for I = ysize(1):-1:1
    for J = ysize(2):-1:xsize(2)+1
      remapdriver([],y{I,J},sprintf([name,'{%1.0f,%1.0f}'],I,J),varID);
    end
  end
elseif any(xsize > ysize)
  % Remove
  csize = ysize;
else
  csize = xsize;
end

% ReMap Common
for I = csize(1):-1:1
  for J = csize(2):-1:1
    remapdriver(x{I,J},y{I,J},sprintf([name,'{%1.0f,%1.0f}'],I,J),varID);
  end
end
end