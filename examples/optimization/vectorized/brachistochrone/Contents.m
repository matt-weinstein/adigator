% ADiGator brachistochrone example w/ fmincon
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% basic_cons.m           - constraint function for use with non-vect
% basic_conswrap.m       - constraint wrapper for use with non-vect
% basic_laggrad.m        - lagrangian gradient file for use with non-vect
% basic_laggradwrap.m    - lagrangian gradient wrapper file for non-vect
% basic_objwrap.m        - objective wrapper for non-vect
% dynamics.m             - dynamics function
% main_basic_1stderivs.m - solve with 1st derivatives, non-vectorized
% main_basic_2ndderivs.m - solve with 2nd derivatives, non-vectorized
% main_noderivs.m        - solve with no supplied derivatives
% main_vect_1stderivs.m  - solve with 1st derivatives, vectorized
% main_vect_2ndderivs.m  - solve with 2nd derivatives, vectorized
% setupproblem.m         - problem setup - initial guess and collocation
% vect_cons.m            - constraint function for use with vect
% vect_conswrap.m        - constraint wrapper for use with vect
% vect_laggrad.m         - lagrangian gradient file for use with vect
% vect_laggradwrap.m     - lagrangian gradient wrapper file for vect
