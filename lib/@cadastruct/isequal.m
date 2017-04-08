function z = isequal(x,y,varargin)
% CADASTRUCT overloaded version of function ISEQUAL (arrays)
% Implementation is only for internal ADiGator algorithm isequal calls,
% logic has not been coded to allow users to use isequal on structures/cell
% arrays.

Dbstuff = dbstack; 
if length(Dbstuff)>1
  CallingFile = Dbstuff(2).file;
else
  CallingFile = [];
end
if (length(CallingFile) > 16 && strcmp(CallingFile(1:16),'adigatortempfunc'))
  error('not yet coded')
elseif nargin > 2
  error('not yet coded')
end

z = isequal(x.val,y.val);