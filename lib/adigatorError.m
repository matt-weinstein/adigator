function adigatorError(IfCount,BroCount,ErrorMsg)
% All user error functions are replaced with this function in the
% intermediate program.
%
% Inputs:
%   IfCount  - integer identifying the conditional if/elseif/else set to
%              which the error belongs
%   BroCount - integer identifying the branch of the conditional set to
%              which the error belongs
%   ErrorMsg - string containing the error message to be printed (won't
%              work well if this is more complex than a simple string)
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR

if ~ADIGATOR.RUNFLAG && ~ADIGATOR.OPTIONS.PREALLOCATE
  ADIGATOR.ERRORLOCS(end+1,:) = [IfCount, BroCount];
elseif ADIGATOR.RUNFLAG == 2 && ~ADIGATOR.EMPTYFLAG
  fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,ErrorMsg,'\n']);
end
end