% ADiGator overloaded @cada class - all doubles of original program
% replaced by objects of cada class in intermediate program
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% abs.m                        - overloaded abs
% acos.m                       - overloaded acos
% acosd.m                      - overloaded acosd
% acosh.m                      - overloaded acosh
% acot.m                       - overloaded acot
% acotd.m                      - overloaded acotd
% acoth.m                      - overloaded acoth
% acsc.m                       - overloaded acsc
% acscd.m                      - overloaded acscd
% acsch.m                      - overloaded acsch
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
% all.m                        - overloaded all
% and.m                        - overloaded and
% any.m                        - overloaded any
% asec.m                       - overloaded asec
% asecd.m                      - overloaded asecd
% asech.m                      - overloaded asech
% asin.m                       - overloaded asin
% asind.m                      - overloaded asind
% asinh.m                      - overloaded asinh
% atan.m                       - overloaded atan
% atan2.m                      - overloaded atan2
% atand.m                      - overloaded atand
% atanh.m                      - overloaded atanh
% cada.m                       - cada constructor function
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
% ceil.m                       - overloaded ceil
% cell.m                       - overloaded cell
% colon.m                      - overloaded colon
% cos.m                        - overloaded cos
% cosd.m                       - overloaded cosd
% cosh.m                       - overloaded cosh
% cot.m                        - overloaded cot
% cotd.m                       - overloaded cotd
% coth.m                       - overloaded coth
% cross.m                      - overloaded cross
% csc.m                        - overloaded csc
% cscd.m                       - overloaded cscd
% csch.m                       - overloaded csch
% ctranspose.m                 - overloaded ctranspose
% diag.m                       - overloaded diag
% display.m                    - overloaded display
% dot.m                        - overloaded dot
% end.m                        - overloaded end
% eq.m                         - overloaded eq
% erf.m                        - overloaded erf
% exp.m                        - overloaded exp
% eye.m                        - overloaded eye
% false.m                      - overloaded false
% fix.m                        - overloaded fix
% floor.m                      - overloaded floor
% full.m                       - overloaded full
% ge.m                         - overloaded ge
% gt.m                         - overloaded gt
% horzcat.m                    - overloaded horzcat
% interp1.m                    - overloaded interp1
% interp2.m                    - overloaded interp2
% interp2old.m                 - overloaded interp2old
% inv.m                        - overloaded inv
% isempty.m                    - overloaded isempty
% isequal.m                    - overloaded isequal
% isequalwithequalnans.m       - overloaded isequalwithequalnans
% ldivide.m                    - overloaded ldivide
% le.m                         - overloaded le
% length.m                     - overloaded length
% log.m                        - overloaded log
% log10.m                      - overloaded log10
% logical.m                    - overloaded logical
% lt.m                         - overloaded lt
% max.m                        - overloaded max
% min.m                        - overloaded min
% minus.m                      - overloaded minus
% mldivide.m                   - overloaded mldivide
% mod.m                        - overloaded mod
% mpower.m                     - overloaded mpower
% mrdivide.m                   - overloaded mrdivide
% mtimes.m                     - overloaded mtimes
% nan.m                        - overloaded nan
% ne.m                         - overloaded ne
% nnz.m                        - overloaded nnz
% nonzeros.m                   - overloaded nonzeros
% not.m                        - overloaded not
% num2str.m                    - overloaded num2str
% numel.m                      - overloaded numel
% ones.m                       - overloaded ones
% or.m                         - overloaded or
% plus.m                       - overloaded plus
% power.m                      - overloaded power
% ppval.m                      - overloaded ppval
% private                      - overloaded priva
% prod.m                       - overloaded prod
% rdivide.m                    - overloaded rdivide
% real.m                       - overloaded real
% rem.m                        - overloaded rem
% repmat.m                     - overloaded repmat
% reshape.m                    - overloaded reshape
% round.m                      - overloaded round
% sec.m                        - overloaded sec
% secd.m                       - overloaded secd
% sech.m                       - overloaded sech
% sign.m                       - overloaded sign
% sin.m                        - overloaded sin
% sind.m                       - overloaded sind
% sinh.m                       - overloaded sinh
% size.m                       - overloaded size
% sparse.m                     - overloaded sparse
% sqrt.m                       - overloaded sqrt
% sub2ind.m                    - overloaded sub2ind
% subsasgn.m                   - overloaded subsasgn
% subsasgnlogical.m            - overloaded subsasgnlogical
% subsasgnold.m                - overloaded subsasgnold
% subsindex.m                  - overloaded subsindex
% subsref.m                    - overloaded subsref
% sum.m                        - overloaded sum
% sum2.m                       - overloaded sum2
% tan.m                        - overloaded tan
% tand.m                       - overloaded tand
% tanh.m                       - overloaded tanh
% times.m                      - overloaded times
% transpose.m                  - overloaded transpose
% true.m                       - overloaded true
% uminus.m                     - overloaded uminus
% uplus.m                      - overloaded uplus
% vertcat.m                    - overloaded vertcat
% xor.m                        - overloaded xor
% zeros.m                      - overloaded zeros
%
% ----------------------------------------------------------------------- %
% DIRECTORIES:
% private                      - private library used by overloaded
%                                procedures