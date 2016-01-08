function y = struct(varargin)
% CADASTRUCT overloaded version of function STRUCT

inputs = varargin;
for i = 2:2:nargin
  if isa(inputs{i},'cadastruct')
    inputs{i} = inputs{i}.val;
  end
end
y = struct(inputs{:});