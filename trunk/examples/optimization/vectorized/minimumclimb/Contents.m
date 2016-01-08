% ADiGator minimum time to climb example w/ fmincon
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% basic_objwrap.m          - objective wrapper
% conswrap.m               - constraint wrapper
% dynamics.m               - dynamics file
% laghesswrap.m            - lagrangian hessian wrapper
% main_1stderivs_nonvect.m - solve with 1st derivatives, non-vect
% main_1stderivs_vect.m    - solve with 1st derivatives, vect
% main_2ndderivs_nonvect.m - solve with 1st derivatives, non-vect
% main_2ndderivs_vect.m    - solve with 2nd derivatives, vect
% setupproblem.m           - get initial guess, set up hermite-simpson 
%                            collocation
