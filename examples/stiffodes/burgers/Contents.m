% ADiGator burgers ODE example
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% burgersfun.m         - original burgers function taken from burgersode.m 
%                        file
% burgersfun_noloop.m  - modified burgers function, (loops removed)
% main.m               - solves burgers ODE supplying derivatives with
%                        ADiGator - ADiGator differentiates the noloop file