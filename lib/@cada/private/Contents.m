% ADiGator overloaded @cada class private library folder
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% cadaCancelDerivs.m    - get rid of a variables derivatives and any
%                         intermediate variable's derivs which are used to
%                         create them
% cadaEmptyEval.m       - empty evaluation, simply propagate .id field and
%                         ADIGATOR.VARINFO.LASTOCC
% cadaRemoveRowsCols.m  - remove rows/cols of an object when have a
%                         dimension mismatch in printing evaluation of loop
% cadaRepDers.m         - repmat derivative elements of a scalar
% cadabinaryarraymath.m - overloaded procedure for all binary array
%                         math operations (e.g. times, plus, minus)
% cadabinarylogical.m   - overloaded procedure for all binary logical array
%                         operations (e.g. ==, ~=)
% cadacreatearray.m     - overloaded procedure for all array instantiation
%                         operations (e.g. zeros, ones)
% cadadername.m         - get derivative variable name for printing
% cadafuncname.m        - get function variable name for printing
% cadaindprint.m        - store indices/generate variable name
% cadainversederiv.m    - print derivative procedures for matrix inverse
% cadamatprint.m        - store data/generate variable name
% cadamtimesderiv.m     - derivative procedures for matrix multiplication
% cadamtimesderivvec.m  - derivative procedures for vectorized matrix
%                         multiplication
% cadaunarylogical.m    - overloaded procedure for all unary logical
%                         operations (e.g. not)
% cadaunarymath.m       - overloaded procedure for all unary array math
%                         operations (e.g. sin, sqrt)
% cadaunion.m           - union of two derivative sparsity patterns
