% ADiGator (Automatic DIfferentiation by GATORs) - A source transformation
% via operator overloading tool for automatic differentiation of MATLAB
% functions.
%
% Copyright 2011-2015 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% website: http://sourceforge.net/projects/adigator/
%
% ----------------------------------------------------------------------- %
% FILES:
% adigator.m                 - main adigator function
% adigatorCreateAuxInput.m   - function for identifying auxiliary numerical
%                              inputs
% adigatorCreateDerivInput.m - function for identifying derivative inputs
% adigatorOptions.m          - function for setting adigator options
%                              structure
% startupadigator.m          - path setup for adigator toolbox
% ----------------------------------------------------------------------- %
% DIRECTORIES:
% doc      - user's guide and licensing information
% examples - various example problems
% lib      - ADiGator library of source transformation routines, overloaded
%            classes, etc
% util     - user utility functions to invoke source transformation