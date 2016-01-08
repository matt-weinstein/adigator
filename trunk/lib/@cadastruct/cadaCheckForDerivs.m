function flag = cadaCheckForDerivs(x)
% This just checks and overloaded structure array to see if it has any
% derivatives
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
flag = structparse(x.val);
end

function flag = structparse(x)

if isstruct(x)
  [M,N] = size(x);
  Fnames = fieldnames(x);
  for I=1:M
    for J = 1:N
      for K = 1:length(Fnames)
        Fstr = Fnames{K};
        flag = structparse(x(I,J).(Fstr));
        if flag
          return
        end
      end
    end
  end
elseif iscell(x)
  [M,N] = size(x);
  for I=1:M
    for J = 1:N
      flag = structparse(x{I,J});
      if flag
        return
      end
    end
  end
elseif isa(x,'cada')
  flag = cadaCheckForDerivs(x);
else
  flag = 0;
end

end