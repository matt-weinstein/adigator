% ADiGator unconstrained minimization example w/ IPOPT
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% dgl2fg.f   - original fortran code computing function and gradient of GL2
%              problem
% gl2f.m     - GL2 objective function file
% gl2main.m  - solves GL2 problem using ADiGator derivatives and IPOPT
% gl2st.m    - initialization of GL2 problem