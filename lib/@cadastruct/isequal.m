function z = isequal(x,y,varargin)
% CADASTRUCT overloaded version of function ISEQUAL (arrays)
Dbstuff = dbstack; CallingFile = Dbstuff(2).file;
if (length(CallingFile) > 16 && strcmp(CallingFile(1:16),'adigatortempfunc'))
  error('not yet coded')
elseif nargin > 2
  error('not yet coded')
end

z = isequal(x.val,y.val);