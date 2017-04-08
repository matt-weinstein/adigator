function varargout = size(x,varargin)
% CADASTRUCT overloaded SIZE routine
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
% just make a dummy CADA variable and use CADA length
NUMvod = ADIGATOR.NVAROFDIFF;
func  = struct('name',x.name,'size',size(x.val),'zerolocs',[],'value',[]);
deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
xdummy = cada(x.id,func,deriv);
if nargout == 1 && nargin == 1
  % This is the case that Matlab Workspace will call - If you put a
  % keyboard in this section you are going to crash matlab if you have the
  % workspace open.
  Dbstuff = dbstack; 
  if length(Dbstuff)>1
    CallingFile = Dbstuff(2).file;
  else
    CallingFile = [];
  end
  if ADIGATOR.OPTIONS.KEYBOARD
    error(['Cannot use size with nargin = nargout = 1 when you have a',...
      '''keyboard'' within the code - causes an error with the MATLAB WorkSpace.']);
  elseif (length(CallingFile) > 16 && strcmp(CallingFile(1:16),'adigatortempfunc'))
    varargout{1} = size(xdummy);
  else
    varargout{1} = func.size;
  end
elseif nargin == 2 && (nargout == 1 || nargout == 0)
  varargout{1} = size(xdummy,varargin{1});
elseif nargin == 1 && (nargout == 2 || nargout == 0)
  [varargout{1}, varargout{2}] = size(xdummy);
elseif nargout > 2
  error('Too many output arguments')
elseif nargin > 2
  error('Too many input arguments')
end