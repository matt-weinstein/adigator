function x = adigatorStructAnalyzer(x,xStr,~)
% function x = adigatorStructAnalyzer(x,xStr,subsflag)
% Non-overloaded version - only called for structures or strings
% 
% Inputs:
% x - variable which was just created
% xStr - string name of variable in the program.
%
% Outputs:
% x - same variable, if structure, contains ordered fields and each field
% is sent to adigatorStructAnalyzer

if isstruct(x) && ~isempty(x)
  x = orderfields(x);
  fnames = fieldnames(x);
  for I = 1:length(fnames)
    F = fnames{I};
    x.(F) = adigatorStructAnalyzer(x.(F),[xStr,'.',F],0);
  end
end

end