% ADiGator user utilities
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% adigator.m                  - main source transformation AD driver
% adigatorColor.m             - CPR coloring scheme
% adigatorCreateAuxInput.m    - creates ADiGator input vars for auxiliary
%                               inputs
% adigatorCreateDerivInput.m  - creates ADiGator input vars for derivative
%                               inputs
% adigatorGenFiles4Fmincon.m  - black box derivative file generation for 
%                               use with fmincon (constrained minimization)
% adigatorGenFiles4Fminunc.m  - black box derivative file generation for 
%                               use with fminunc 
%                               (unconstrained minimization)
% adigatorGenFiles4Fsolve.m   - black box derivative file generation for 
%                               use with fsolve (system of equations)
% adigatorGenFiles4gpops2.m   - black box derivative file generation for 
%                               use with GPOPS2 optimal control software
% adigatorGenFiles4Ipopt.m    - black box derivative file generation for 
%                               use with IPOPT (constrained minimization)
% adigatorGenHesFile.m        - generates first and second derivative files 
%                               with wrapper files
% adigatorGenJacFile.m        - generates first derivative files with 
%                               wrapper file
% adigatorOptions.m           - creates adigator Options structure
% adigatorProjectVectLocs.m   - projects locations of non-vectorized 
%                               Jacobian into those of vectorized Jacobian
% adigatorUncompressJac.m     - uncompresses Jacobians for use with CPR
%                               coloring
% ----------------------------------------------------------------------- %
% DIRECTORIES:
% private                     - adigatortempfunc#.m files get stored here                    