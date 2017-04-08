function x = adigatorCreateAuxInput(xsize,value)
% ADiGator auxiliary input variable creation routine.
%
% --------------------------- Usage ------------------------------------- %
%  x = adigatorCreateAuxInput(xsize)
% This function creates an input,x, (which has NO derivatives, but is not 
% a fixed value), that is to be used with the function source-to-derivative
% source transformation function, adigator.
%             OR
%  x = adigatorCreateAuxInput(xsize,value)
% This is useful if you have a vectorized input which is actually a fixed
% value. If the first dimension is vectorized, then value should be a row
% vector, if second dimension is vectorized, then value should be a column
% vector.
%
% ------------------------ Input Information ---------------------------- %
% xsize - size of auxiliary input variable
% value - fixed value of auxiliary input variable (optional, should ony be
% used with vectorized, known auxiliary inputs)
% 
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% See also: adigatorCreateDerivInput adigatorOptions adigator

if isequal(size(xsize),[1 2]) && isequal(xsize,floor(xsize))
  func.size = xsize;
else
  error('first input xsize must be integer array of size 1 by 2')
end
if isinf(xsize(1)) && isinf(xsize(2))
  error('only one dimension of the input may be vectorized');
end

if nargin == 2
  valsize = size(value);
  if any(xsize(~isinf(xsize)) ~= valsize(~isinf(xsize))) || any(valsize(isinf(xsize))~=1)
    error('Invalue value input');
  end
  func.value = value;
end

x = adigatorInput(func,[]);
