function x = adigatorStructAnalyzer(x,xStr,~)
% function x = adigatorStructAnalyzer(x,xStr,subsflag)
% Non-overloaded version - only called for structures or strings

if isstruct(x) && ~isempty(x)
  x = orderfields(x);
  fnames = fieldnames(x);
  for I = 1:length(fnames)
    F = fnames{I};
    x.(F) = adigatorStructAnalyzer(x.(F),[xStr,'.',F],0);
  end
end

end