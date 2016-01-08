function z = cadaUnionVars(x,y)
% Union two cadastruct vars together
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

if ~isa(x,'cadastruct')
  z = y;
  xval = x;
  yval = y.val;
elseif ~isa(y,'cadastruct')
  xval = x.val;
  yval = y;
  z = x;
else
  xval = x.val;
  yval = y.val;
  z = x;
end

z.val = uniondriver(xval,yval,z.name);
end

function z = uniondriver(x,y,name)
if isa(x,'cadastruct')
  x = x.val;
end
if isa(y,'cadastruct')
  y = y.val;
end
xclass = class(x); yclass = class(y);
if ~strcmp(xclass,yclass)
  if (isa(x,'cada') && prod(x.func.size)  == 0) || (isnumeric(x) && isempty(x))
    z = y;
  elseif isa(y,'cada') && prod(y.func.size) == 0 || (isnumeric(y) && isempty(y))
    z = x;
  else
    error([name,' takes on either an ',xclass,' or a ',yclass,...
      ' depending upon a loop iteration, conditional statement,',...
      ' function call, etc. This is not allowed']);
  end
else
  switch xclass
    case 'cada'
      z = cadaUnionVars(x,y);
    case 'struct'
      z = structureunion(x,y,name);
    case 'cell'
      z = cellunion(x,y,name);
    otherwise
     z = x; 
  end
end

end

function z = structureunion(x,y,name)
xMrow = size(x,1); xNcol = size(x,2); xfields = fieldnames(x); 
yMrow = size(y,1); yNcol = size(y,2); yfields = fieldnames(y);
zMrow = max(xMrow,yMrow); zNcol = max(xNcol,yNcol);

% Need to determine 4 sets {i,j,f}, where z(i,j).f = 
% Set 1: x(i,j).f U y(i,j).f
% Set 2: x(i,j).f;
% Set 3: y(i,j).f;
% Set 4: empty
zfields = unique([xfields;yfields]);
zset = [zfields.'; repmat({cell(zMrow,zNcol)},1,length(zfields))];
% Initialize z
z    = struct(zset{:});
if zMrow*zNcol == 1
  for Fcount = 1:length(zfields)
    F = zfields{Fcount};
    xfieldflag = any(strcmp(F,xfields));
    yfieldflag = any(strcmp(F,yfields));
    if xfieldflag
      if yfieldflag
        % Case 1
        z.(F) = uniondriver(x.(F),y.(F),sprintf('%s.%s',name,F));
      else
        % Case 2
        z.(F) = x.(F);
      end
    elseif yfieldflag
      % Case 3
      z.(F) = y.(F);
    end
  end
else
  for Fcount = 1:length(zfields)
    F = zfields{Fcount};
    xfieldflag = any(strcmp(F,xfields));
    yfieldflag = any(strcmp(F,yfields));
    for I = 1:zMrow
      for J = 1:zNcol
        if xfieldflag && I <= xMrow && J <= xNcol
          if yfieldflag && I <= yMrow && J <= yNcol
            % Case 1
            z(I,J).(F) = uniondriver(x(I,J).(F),y(I,J).(F),sprintf('%s(%1.0f,%1.0f).%s',name,I,J,F));
          else
            % Case 2
            z(I,J).(F) = x(I,J).(F);
          end
        elseif yfieldflag && I <= yMrow && J <= yNcol
          % Case 3
          z(I,J).(F) = y(I,J).(F);
        end
      end
    end
  end
end
end

function z = cellunion(x,y,name)
xMrow = size(x,1); xNcol = size(x,2);
yMrow = size(y,1); yNcol = size(y,2);
zMrow = max(xMrow,yMrow); zNcol = max(xNcol,yNcol);

z = cell(zMrow,zNcol);
for I = 1:zMrow
  for J = 1:zNcol
    if I <= xMrow && J <= xNcol
      if I <= yMrow && J <= yNcol
        z{I,J} = uniondriver(x{I,J},y{I,J},sprintf('%s{%1.0f,%1.0f}',name,I,J));
      else
        z{I,J} = x{I,J};
      end
    elseif I <= yMrow && J <= yNcol
      z{I,J} = y{I,J};
    end
  end
end
end