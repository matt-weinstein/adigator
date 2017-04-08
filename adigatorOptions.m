function options = adigatorOptions(varargin)
% ADiGator Options Structure Generator: Called in order set/modify default
% ADiGator options.
% 
% ------------------------------ Usage ---------------------------------- %
% options = adigatorOptions(field1,value1,field2,value2,...)
%
% ------------------------ Input Information ---------------------------- %
% field:  string name of option to be modified
% value:  value to set option to
%
% ----------------------- Output Information ---------------------------- %
% options: structure to be passed to adigator, or any of the utility
%          ADiGator generation routines (adigatorGenJacFile,
%          adigatorGenHesFile, adigatorGenFiles4Fmincon,
%          adigatorGenFiles4Ipopt, adigatorGenFiles4gpops2)
%
%                                OPTIONS:  
% ------------------------------------------------------------------------
%    AUXDATA:  1 - auxiliary inputs will always have the same sparsity
%                  pattern but may change numeric values 
%              0 - auxiliary inputs will always have the same numeric
%                  values (required if inputs are a size to be looped on or
%                  a reference index, etc.) (default)
%       ECHO:  1 - echo to screen the transformation progress (default)
%              0 - dont echo
%     UNROLL:  1 - unroll loops and sub functions in derivative program
%              0 - keep loops and sub functions rolled in derivative 
%                  (default)
%   COMMENTS:  1 - print comments to derivative file giving the lines of
%                  user code which correspond to the printed statements
%                  (default)
%              0 - dont print comments
%  OVERWRITE:  1 - if the user supplies a derivative file name which
%                  corresponds to an already existing file, then setting
%                  this option will overwrite the existing file. (default
%                  for adigatorGenFiles4gpops2, adigatorGenJacFile,
%                  adigatorGenHesFile, adigatorGenFiles4Ipopt,
%                  adigatorGenFiles4Fmincon)
%              0 - if a file already exists with the same name as the given
%                  derivative file name, then it will not overwrite and
%                  error out instead (default for adigator)
% MAXWHILEITER:k - maximum number of iterations to attempt to find a static
%                  input/output for WHILE loops (default set to 10) -
%                  increasing this will increase derivative file generation
%                  times when using WHILE loops
% COMPLEX:     0 - do not expect any variables to be complex, use 
%                  non-complex forms of abs, ctranspose, dot (default)
%              1 - expect variables to be complex, use complex forms of
%                  ctranspose, abs, dot.
% ------------------------------------------------------------------------
%
% NOTES:    The default value of the OVERWRITE option changes depending
%           upon whether the basic adigator file is being called or one of
%           the wrapper generation files.
%
% If desired, the defaults of each option may be changed by editing this
% file.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% See also adigator, adigatorGenFiles4gpops2

% Set Defaults
options.auxdata   = 0;
options.echo      = 1;
options.unroll    = 0;
options.comments  = 1;
options.overwrite = 0;
options.optoutput = 0;
options.maxwhileiter = 10;
options.complex   = 0;
if nargin/2 ~= floor(nargin/2)
  error('Inputs to adigatorOptions must come in field/value pairs')
end

% Set user wanted options
for i = 1:nargin/2
  field = lower(varargin{2*(i-1)+1});
  value = varargin{2*i};
  switch field
    case {'auxdata','echo','unroll','comments','overwrite','genpat',...
        'optoutput','complex'}
      options.(field) = logical(value);
    case 'maxwhileiter'
      options.(field) = value;
    otherwise
      warning(['Invalid option field: ',field])
  end
end