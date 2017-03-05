function adigatorIfLooperi(ifLoopi,IfCount)
% This is a transformation routine that is used when we are unrolling loops
% but encounter an IF statement. It is placed within the "if loop" set by
% adigatorIfLooper.
%
% Inputs: 
%   ifLoopi  - variable set by adigatorIfLooper which determines what type
%              evaluation of the conditional statement must be performed.
%   IfCount  - integer identifying the conditional if/elseif/else set
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.OPTIONS.PREALLOCATE
  return
end
ADIGATOR.RUNFLAG = ifLoopi;
if ifLoopi
  ADIGATOR.VARINFO.COUNT = ADIGATOR.IFDATA(IfCount).BROS(1).START;
  ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  if ifLoopi == 1
    ADIGATOR.PRINT.FLAG = 0;
  else
    ADIGATOR.PRINT.FLAG = 1;
  end
end