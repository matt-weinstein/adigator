function output = adigatorGenJacFile(UserFunName,UserFunInputs,varargin)
% ADiGator Jacobian File Generation Function: this function is used when
% you wish to generate a Jacobian of a function with an input variable of
% differentiation (and any auxiliary inputs), and a single output. This
% simply calls the function adigator and creates a wrapper function s.t.
% the input of the resulting file is the same as the input to the original
% user function, but outputs both the Jacobian and the original output.
%
% ------------------------------ Usage -----------------------------------
% output = adigatorGenJacFile(UserFunName,UserFunInputs)
%                   or
% output = adigatorGenJacFile(UserFunName,UserFunInputs,Options)
%
% ------------------------ Input Information -----------------------------
% UserFunName: String name of the user function to be differentiated
%
% UserFunInputs: N x 1 cell array containing the inputs to the UserFun
%                - the input (or cell array element/structure field)
%                corresponding to the variable of differentiation must be
%                created using the adigatorCreateDerivInput function.
%                i.e. if the first input is the variable of
%                differentiation, then the first input must be created
%                using adigatorCreateDerivInput prior to calling adigator.
%                - any other numeric inputs should be defined as they will
%                be when calling the derivative function. These will be
%                assumed to have fixed sizes and zero locations, but the
%                non-zero locations may change values. If the values are
%                always fixed, then adigatorOptions may be used to change
%                the handling of these auxiliary inputs.
%                - auxiliary inputs may also be created using the
%                adigatorCreateAuxInput function.
%
% Options (optional): option structure generated by adigatorOptions
%                     function
%
% ------------------------ Output Information ----------------------------
% Assuming UserFunName = 'myfun', then the Jacobian is named 'myfun_Jac'. 
% The output of adigatorGenJacFile is the structure:
%     output.FunctionFile = 'myfun'
%     output.JacobianFile = 'myfun_Jac'
%     output.JacobianStructure = sparse ones and zeros.
% The generated Jacobian file has the same input structure as the original
% user function. The output of Jacobian file is [Jac, Fun].
%
% ----------------------- Additional Information -------------------------
% The Jacobian is built as a sparse matrix under the condition that
% numel(Jac) >= 250 & nnz(Jac)/numel(Jac) <= 3/4, otherwise it is built as
% a full matrix.
%
% If the output is y, input is x, numel(y) = m, numel(x) = n...
%     case: n > 1, size(Jac) = [m n]
%     case: n = 1, size(Jac) = size(y)
%
% The functions generated Jacobian file is simply wrapper files for the
% ADiGator generated file (named 'myfun_ADiGatorJac')
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% see also adigator, adigatorCreateDerivInput, adigatorCreateAuxInput,
% adigatorOptions, adigatorGenHesFile

if ~ischar(UserFunName)
  error(['First input to adigator must be string name of function to be ',...
    'differentiated']);
end
JacFileName    = [UserFunName,'_Jac'];         % Name of wrapper file
AdiJacFileName = [UserFunName,'_ADiGatorJac']; % Name of derivative file

% Options
if nargin == 2
  opts.overwrite = 1;
else
  opts = varargin{1};
  if ~isfield(opts,'overwrite')
    opts.overwrite = 1;
  end
end

% Quick input check
if ~iscell(UserFunInputs)
  error(['Second input to adigator must be cell array of inputs to ',...
    'the function described by first input string']);
end

% Find derivative input
derflag = 0;
for I = 1:numel(UserFunInputs)
  x = UserFunInputs{I};
  if isa(x,'adigatorInput')
    if ~isempty(x.deriv)
      if derflag > 0
        error('adigatorGenJacFile is only used for single derivative variable input')
      end
      derflag = I;
    end
    if any(isinf(x.func.size))
      error('adigatorGenJacFile not written for vectorized functions')
    end
  end
end
if derflag == 0
  error('derivative input of user function not found - possibly embedded within a cell/structure, use adigator function if this is the case');
end
UserFun = str2func(UserFunName);
% Output Check
if nargout(UserFun) ~= 1
  error('User function must contain single output');
end


% File checks
CallingDir = cd;
if exist([CallingDir,filesep,JacFileName,'.m'],'file');
  if opts.overwrite
    delete([CallingDir,filesep,JacFileName,'.m']);
    rehash
  else
    error(['The file ',CallingDir,filesep,JacFileName,'.m already exists, ',...
      'quitting transformation. To set manual overwrite of file use ',...
      '''''adigatorOptions(''OVERWRITE'',1);''''. Alternatively, delete the ',...
      'existing file and any associated .mat file.']);
  end
end


% Call adigator
[adiout,FunctionInfo] = adigator(UserFunName,UserFunInputs,AdiJacFileName,opts);
adiout = adiout{1};
fid = fopen([JacFileName,'.m'],'w+');

InputStrs = FunctionInfo.Input.Names.';
xstr = InputStrs{derflag};
for I = 1:length(InputStrs)
  InputStrs{I} = [InputStrs{I},','];
end
InputStr1 = cell2mat(InputStrs);
InputStr1(end) = [];

functionstr = ['function [Jac,Fun] = ',JacFileName,'(',InputStr1,')\n'];

fprintf(fid,['%% ',functionstr,]);
fprintf(fid,'%% \n');

% Print Function Header
fprintf(fid,'%% Jacobian wrapper file generated by ADiGator\n');
fprintf(fid,['%% ',char(169),'2010-2014 Matthew J. Weinstein and Anil V. Rao\n']);
fprintf(fid,'%% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ \n');
fprintf(fid,'%% Contact: mweinstein@ufl.edu\n');
fprintf(fid,'%% Bugs/suggestions may be reported to the sourceforge forums\n');
fprintf(fid,'%%                    DISCLAIMER\n');
fprintf(fid,'%% ADiGator is a general-purpose software distributed under the GNU General\n');
fprintf(fid,'%% Public License version 3.0. While the software is distributed with the\n');
fprintf(fid,'%% hope that it will be useful, both the software and generated code are\n');
fprintf(fid,'%% provided ''AS IS'' with NO WARRANTIES OF ANY KIND and no merchantability\n');
fprintf(fid,'%% or fitness for any purpose or application.\n\n');

fprintf(fid,functionstr);
% Change the derivative input..
x = UserFunInputs{derflag};
vodname = x.deriv.vodname;
xfunstr = ['gator_',xstr,'.f'];
xderstr = ['gator_',xstr,'.d',vodname];
fprintf(fid,[xfunstr,' = ',xstr,';\n']);
fprintf(fid,[xderstr,' = ones(%1.0f,1);\n'],prod(x.func.size));

InputStrs{derflag} = ['gator_',xstr,','];
InputStr2 = cell2mat(InputStrs);
InputStr2(end) = [];

ystr = FunctionInfo.Output.Names{1};
% Call the ADiGatorJac file
fprintf(fid,[ystr,' = ',AdiJacFileName,'(',InputStr2,');\n']);

% Check to see how many non-zeros in Jacobian
xsize = x.func.size;
ysize = adiout.func.size;
dydxsize = [prod(ysize), prod(xsize)];
dydxnumel  = dydxsize(1)*dydxsize(2);
if dydxsize(1) == 1 && all(xsize>1)
  dydxsize = xsize;
  ysize = [xsize(1) 1];
  xsize = [xsize(2) 1];
elseif dydxsize(2) == 1 && all(ysize>1)
  dydxsize = ysize;
  xsize = [ysize(2) 1];
  ysize = [ysize(1) 1];
end
dydxnnz  = size(adiout.deriv.nzlocs,1);
% If dydx has => 250 elements and has <= 75% nonzeros, project into sparse
% matrix, otherwise project into full matrix.
dy = [ystr,'.d',vodname];
if dydxnnz == dydxnumel
  fprintf(fid,['Jac = reshape(',dy,',[%1.0f %1.0f]);\n'],dydxsize);
elseif dydxsize(1) == 1 && dydxsize(2) == 1
  fprintf(fid,['Jac = ',dy,';\n']);
elseif dydxsize(1) == 1
  fprintf(fid,'Jac = zeros(1,%1.0f);',dydxsize(2));
  fprintf(fid,['Jac(',dy,'_location) = ',dy,';\n']);
elseif dydxsize(2) == 1
  fprintf(fid,'Jac = zeros(%1.0f,1);',dydxsize(2));
  fprintf(fid,['Jac(',dy,'_location) = ',dy,';\n']);
else
  dyloc = [dy,'_location'];
  if ~any(ysize == 1)
    % Output is matrix
    fprintf(fid,['funloc = (',dyloc,'(:,2)-1)*%1.0f + ',dyloc,'(:,1);\n'],ysize(1));
    rowstr = 'funloc';
    if ~any(xsize == 1)
      % Input is matrix
      fprintf(fid,['varloc = (',dyloc,'(:,4)-1)*%1.0f + ',dyloc,'(:,3);\n'],xsize(1));
      colstr = 'varloc';
    else
      colstr = [dyloc,'(:,3)'];
    end
  else
    rowstr = [dyloc,'(:,1)'];
    if ~any(xsize == 1)
      % Input is matrix
      fprintf(fid,['varloc = (',dyloc,'(:,3)-1)*%1.0f + ',dyloc,'(:,2);\n'],xsize(1));
      colstr = 'varloc';
    else
      colstr = [dyloc,'(:,2)'];
    end
  end
  if dydxnumel >= 250 && dydxnnz/dydxnumel <= 3/4
    % Project Sparse
    fprintf(fid,['Jac = sparse(',rowstr,',',colstr,',',dy,',%1.0f,%1.0f);\n'],dydxsize);
  else
    % Project Full
    fprintf(fid,'Jac = zeros(%1.0f,%1.0f);\n',dydxsize);
    fprintf(fid,['Jac((',colstr,'-1)*%1.0f+',rowstr,') = ',dy,';\n'],dydxsize(1));
  end
end
fprintf(fid,['Fun = ',ystr,'.f;\n']);
fprintf(fid,'end');
fclose(fid);
rehash

output.FunctionFile = UserFunName;
output.JacobianFile = JacFileName;
output.JacobianStructure = sparse(adiout.deriv.nzlocs(:,1),...
  adiout.deriv.nzlocs(:,2),ones(dydxnnz,1),dydxsize(1),dydxsize(2));

fprintf(['\n<strong>adigatorGenJacFile</strong> successfully generated Jacobian wrapper file: ''',JacFileName,''';\n\n']);