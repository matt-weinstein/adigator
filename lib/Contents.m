% ADiGator source transformation library
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% CLASSES:
% cada                          - primary overloaded class, all doubles of
%                                 original program replaced by objects of
%                                 cada class in intermediate program
% cadastruct                    - all cells/structures of original program
%                                 replaced by objects of cadastruct class
%                                 in intermediate program, cadastruct
%                                 objects contain cell/structure arrays
%                                 which also contain cada objects.
% adigatorInput                 - class to define adigator input variables
%
% ----------------------------------------------------------------------- %
% FILES:
% adigatorAssignImpVarNames.m   - adds new row to ADIGATOR.VARINFO.NAMELOCS
%                                 when new user variable encountered
% adigatorAssignOvermapScheme.m - uses data collected in ADIGATOR.VARINFO
%                                 on parsing run to assign overmapping 
%                                 scheme
% adigatorBreakCont.m           - called from intermediate program where a
%                                 break/continue statement used to be     
% adigatorError.m               - called from intermediate program where an
%                                 error statement used to be
% adigatorEvalInterp2pp.m       - evaluates a 2D piecewise polynomial as
%                                 built by adigatorGenInterp2pp
% adigatorForInitialize.m       - initialization of for loops in
%                                 intermediate program
% adigatorForIterEnd.m          - called at end of each loop iteration
% adigatorForIterStart.m        - called at start of each loop iteration
% adigatorFunctionEnd.m         - called at end of each adigatortempfunc#
% adigatorFunctionInitialize.m  - called at start of each adigatortempfunc#
% adigatorGenInterp2pp.m        - generates coefficients of 2D
%                                 interpolation polynomials
% adigatorIfInitialize.m        - initialization of conditional
%                                 if/elseif/else statements in intermediate
%                                 program
% adigatorIfIterEnd.m           - called at start of each if/elseif/else
%                                 branch
% adigatorIfIterStart.m         - called at end of each if/elseif/else
%                                 branch
% adigatorIfLooper.m            - initializes a loop on an if statement if
%                                 ADIGATOR.OPTIONS.UNROLL = 1
% adigatorIfLooperi.m           - called on each iteration of an if 
%                                 statement if ADIGATOR.OPTIONS.UNROLL = 1
% adigatorMakeNumeric.m         - transforms a double object into a known
%                                 numeric cada object
% adigatorPrintTempFiles.m      - transforms user functions into
%                                 intermediate functions
% adigatorSeperateFunLines.m    - looks for multiple lines of code on
%                                 single function line
% adigatorSetCellEvalFlag.m     - sets a global flag so that any celleval
%                                 commands in the intermediate program do
%                                 not invoke cadastruct operations
% adigatorStructAnalyzer.m      - recursively parses cells/structures
%                                 (non-overloaded version)
% adigatorVarAnalyzer.m         - analyzes variables after they have been
%                                 assigned in intermediate program 
%                                 (non-overloaded version)
% ----------------------------------------------------------------------- %
% DIRECTORIES:
% cadaUtils                     - common utility functions called by
%                                 cada/cadastruct routines