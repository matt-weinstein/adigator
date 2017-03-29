% ADiGator DCAL controller ODE example
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% ----------------------------------------------------------------------- %
% FILES:
% main.m       - sets up the DCAL controller problem for 2 link robot -
%                differentiates getqd, getYd, and TwoLinkSys with ADiGator
% getqd.m      - desired trajectory file
% getYd.m      - desired trajectory matrix file
% TwoLinkSys.m - robot dynamics with DCAL controller