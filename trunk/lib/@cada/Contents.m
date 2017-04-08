% ADiGator overloaded @cada class 
% All numeric objects of original program replaced by objects of cada class
% in intermediate program during ADiGator source transformation.
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% adigatorAnalyzeForData.m     - looks at all organizational operation data
%                                collected for FOR loops, gets data
%                                dependencies, determines iteration
%                                dependent organizational operation
%                                derivative routines
% adigatorEvalInterp2pp.m      - overloaded adigatorEvalInterp2pp
% adigatorPrintOutputIndices.m - takes the deriv.nzlocs field and prints
%                                out derivative matrix sizes/indices (so
%                                that users can map non-zero derivatives
%                                into proper matrices)
% adigatorStructAnalyzer.m     - recursively parses cells/structures
%                                (cada overloaded version)
% adigatorVarAnalyzer.m        - analyzes variables after they have been
%                                assigned in intermediate program 
%                                (cada overloaded version)
% cada.m                       - cada classdef file
% cadabinaryarraymath.m        - overloaded procedure for all binary array
%                                math operations (e.g. times, plus, minus)
% cadabinarylogical.m          - overloaded procedure for all binary 
%                                logical array operations (e.g. ==, ~=)
% cadacreatearray.m            - overloaded procedure for all array 
%                                instantiation operations (e.g. zeros)
% cadaunarylogical.m           - overloaded procedure for all unary logical
%                                operations (e.g. not)
% cadaunarymath.m              - overloaded procedure for all unary array 
%                                math operations (e.g. sin, sqrt)
% cadaCheckForDerivs.m         - check to see if object has derivatives
%                                (cada overloaded version)
% cadaOverMap.m                - builds overmapped objects in overmapping
%                                evaluation, stores variables when
%                                necessary, calls cadaPrintReMap when
%                                necessary (cada overloaded version)
% cadaPrintReMap.m             - prints re-mapping procedures when either
%                                going from an undermapped object to an
%                                overmapped object (adding zeros), or the
%                                opposite (removing zeros) (cada overloaded
%                                version)
% cadaUnionVars.m              - cada overloaded union operator 
% colon.m                      - overloaded colon
% cross.m                      - overloaded cross
% diag.m                       - overloaded diag
% horzcat.m                    - overloaded horzcat
% interp1.m                    - overloaded interp1
% interp2.m                    - overloaded interp2
% inv.m                        - overloaded inv
% isempty.m                    - overloaded isempty
% isequal.m                    - overloaded isequal
% isequalwithequalnans.m       - overloaded isequalwithequalnans
% length.m                     - overloaded length
% max.m                        - overloaded max
% min.m                        - overloaded min
% mldivide.m                   - overloaded mldivide
% mpower.m                     - overloaded mpower
% mrdivide.m                   - overloaded mrdivide
% mtimes.m                     - overloaded mtimes
% nnz.m                        - overloaded nnz
% nonzeros.m                   - overloaded nonzeros
% numel.m                      - overloaded numel
% ppval.m                      - overloaded ppval
% prod.m                       - overloaded prod
% repmat.m                     - overloaded repmat
% reshape.m                    - overloaded reshape
% size.m                       - overloaded size
% sparse.m                     - overloaded sparse
% sub2ind.m                    - overloaded sub2ind
% subsasgn.m                   - overloaded subsasgn
% subsref.m                    - overloaded subsref
% sum.m                        - overloaded sum
% transpose.m                  - overloaded transpose
% vertcat.m                    - overloaded vertcat
%
% ----------------------------------------------------------------------- %
% PRIVATE DIRECTORY FILES:
% cadaCancelDerivs.m    - get rid of a variables derivatives and any
%                         intermediate variable's derivs which are used to
%                         create them
% cadainversederiv.m    - print derivative procedures for matrix inverse
% cadamtimesderiv.m     - derivative procedures for matrix multiplication
% cadamtimesderivvec.m  - derivative procedures for vectorized matrix
%                         multiplication
% cadaRemoveRowsCols.m  - remove rows/cols of an object when have a
%                         dimension mismatch in printing evaluation of loop
% cadaRepDers.m         - repmat derivative elements of a scalar
% cadaunion.m           - union of two derivative sparsity patterns