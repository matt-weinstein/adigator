function adigatorSetCellEvalFlag(flagval)
% function adigatorSetCellEvalFlag(flagval)
% This function was created with V2 when structures started getting handled
% differently - it is used to tell @cadastruct/subsasgn when a
% cellfun(@eval,cellstr) is being used to modify the overloaded workspace.
% It simply sets a global flag equal to the input.

global ADIGATOR
ADIGATOR.CELLEVALFLAG = flagval;
end