% ADiGator overloaded @cadastruct class - all cell/structs of original
% program replaced by objects of cadastruct class in intermediate program
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% adigatorPrintOutputIndices.m - print output derivative indices, calls
% overloaded cada version
% adigatorStructAnalyzer.m     - recursively parses cells/structures
%                                (overloaded cadastruct version)
% adigatorVarAnalyzer.m        - analyzes variables after they have been
%                                assigned in intermediate program 
%                                (overloaded cadastruct version)
% cadaCheckForDerivs.m         - check to see if any of the cada objects
%                                belonging to the cell/structure have
%                                derivatives
% cadaGetStruct.m              - extracts the actual structure/cell array
% cadaGetStructID.m            - extracts the id of the object
% cadaIsArray.m                - returns value of .arrayflag field
% cadaOverMap.m                - builds overmapped objects in overmapping
%                                evaluation, stores variables when
%                                necessary, calls cadaPrintReMap when
%                                necessary (cadastruct overloaded version)
% cadaPrintReMap.m             - prints re-mapping procedures when either
%                                going from an undermapped object to an
%                                overmapped object (adding zeros), or the
%                                opposite (removing zeros) (cadastruct
%                                overloaded version)
% cadaUnionVars.m              - cadastruct overloaded union operator 
% cadastruct.m                 - cadastruct instantiator/class defs
% cadastructDecomp.m           - decompose a cadastruct object
% ctranspose.m                 - overloaded ctranspose
% display.m                    - overloaded display
% horzcat.m                    - overloaded horzcat
% isequal.m                    - overloaded isequal
% length.m                     - overloaded length
% numel.m                      - overloaded numel
% ppval.m                      - overloaded ppval
% private                      - overloaded priva
% repmat.m                     - overloaded repmat
% reshape.m                    - overloaded reshape
% size.m                       - overloaded size
% struct.m                     - overloaded struct
% subsasgn.m                   - overloaded subsasgn
% subsref.m                    - overloaded subsref
% transpose.m                  - overloaded transpose
% vertcat.m                    - overloaded vertcat
% ----------------------------------------------------------------------- %
% DIRECTORIES:
% private                      - private library used by overloaded
%                                procedures