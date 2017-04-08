function varargout = adigatorVarAnalyzer(FunString,varargin)
% This module is the NON-Overloaded Version of adigatorVarAnalyzer. 
%
% This is called from the temporary functions after a line of user code has
% been evaluated in order to analyze the outputs. As this is the
% non-overloaded version, it will only be called if all inputs are
% non-overloaded.
% -----------------------Input Information------------------------------- %
% FunString - the actual User's line of code which has just been evaluated
% varargin  - Information on all of the outputs from the user line of code
% which has just been evaluated. For each output, the inputs to
% adigatorVarAnalyzer are:
%         1. actual output
%         2. string of the name of the output (as defined by user)
%         3. flag stating whether output was subsasgn'd or not
% -----------------------Output Information------------------------------ %
% varargout - overloaded outputs, may or may not have some properties
% changed (may also be cells/structures of overloaded objects)
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
% ----------------------------- Parse Inputs ---------------------------- %
NUMvars   = nargout;
varargout = cell(NUMvars,1);
if strcmp(FunString,'global')
  % Global variables have a special case.
  VarStrings = varargin;
  Variables  = cadaGetGlobalVars(VarStrings);
  SubsFlags  = zeros(NUMvars,1);
else
  Variables  = cell(NUMvars,1);
  VarStrings = cell(NUMvars,1);
  SubsFlags  = zeros(NUMvars,1);
  for Vcount = 1:NUMvars
    Variables{Vcount}  = varargin{1+(Vcount-1)*3};
    VarStrings{Vcount} = varargin{2+(Vcount-1)*3};
    SubsFlags(Vcount)  = varargin{3+(Vcount-1)*3};
  end
end
if ADIGATOR.OPTIONS.PREALLOCATE
  varargout = Variables;
  return
end
PreOpCount = ADIGATOR.PREOPCOUNT;

if ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  if NUMvars
    % Is some sort of variable assignment.
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if isnumeric(x)
        % need to give this guy his own Operation
        % Count, and turn it into a symbolic
        x = adigatorMakeNumeric(x);
        ADIGATOR.VARINFO.LASTOCC(ADIGATOR.VARINFO.COUNT,1)...
          = ADIGATOR.VARINFO.COUNT;
        % Set his VarName
        adigatorAssignImpVarNames(ADIGATOR.VARINFO.COUNT,VarStrings{Vcount},SubsFlags(Vcount));
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      elseif isstruct(x) || iscell(x)
        % Is a cell or structure - Need to give a new operation count to
        % each numeric/overloaded object embedded in the cell/structure.
        % Numeric cells/fields need to be turned into Overloaded Objects.
        if PreOpCount < ADIGATOR.VARINFO.COUNT && Vcount == 1
          xID = ADIGATOR.VARINFO.COUNT;
          ADIGATOR.VARINFO.NAMELOCS(PreOpCount:xID-1,2)=1:xID-PreOpCount;
        end
        if ADIGATOR.VARINFO.COUNT > PreOpCount && ~isempty(strfind(FunString,'struct'))
          error(['Cannot parse statement:''',FunString,''' due to overloaded operation ',...
            'being performed within a cell/structure instantiation']);
        end
        x = structParse(x,VarStrings{Vcount},0);
        if isstruct(x)
          if ADIGATOR.VARINFO.COUNT == PreOpCount
            % Give this variable his own count - might need for structure
            % array
            x = cadastruct(x,VarStrings{Vcount},ADIGATOR.VARINFO.COUNT,0);
            adigatorAssignImpVarNames(ADIGATOR.VARINFO.COUNT,VarStrings{Vcount},SubsFlags(Vcount));
            ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
          else
            x = cadastruct(x,VarStrings{Vcount},[],0);
          end
        end
      else
        % String
      end
      varargout{Vcount} = x;
    end
    if strcmp(FunString,'global')
      ADIGATOR.VARINFO.NAMELOCS(PreOpCount...
        :ADIGATOR.VARINFO.COUNT-1,3) = -Inf;
    end
  elseif strcmp(FunString,'keyboard')
    
  elseif strcmp(FunString,'return')
    
  end
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  if NUMvars
    % Is some sort of variable assignment.
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if isnumeric(x)
        % Variable is Numeric
        if all(size(x)) == 0
          % Check to see if need to make a cadastruct array
          xid = ADIGATOR.VARINFO.COUNT;
          xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
          asgnlocs = find(ADIGATOR.VARINFO.NAMELOCS(:,1)==xnameloc);
          nextasgn = asgnlocs(asgnlocs>xid);
          if ~isempty(nextasgn) && ~ADIGATOR.VARINFO.NAMELOCS(nextasgn(1),3)
            x = cadastruct([],VarStrings{Vcount},xid,1);
          else
            x = adigatorMakeNumeric(x);
          end
        else
          x = adigatorMakeNumeric(x);
        end
        % OverMapping
        x = cadaOverMap(x);
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      elseif isstruct(x) || iscell(x)
        % Cell or Structure
        x = structParse(x,VarStrings{Vcount},0);
        if isstruct(x)
          if ADIGATOR.VARINFO.COUNT == PreOpCount
            % This guy gets his own variable count
            % Check to see if need to make a cadastruct array
            xid = ADIGATOR.VARINFO.COUNT;
            xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
            asgnlocs = find(ADIGATOR.VARINFO.NAMELOCS(:,1)==xnameloc);
            nextasgn = asgnlocs(asgnlocs>xid);
            if ~isempty(nextasgn) && ~ADIGATOR.VARINFO.NAMELOCS(nextasgn(1),3)
              x = cadastruct([],VarStrings{Vcount},xid,1);
              x = cadaOverMap(x);
            else
              x = cadastruct([],VarStrings{Vcount},xid,0);
            end
            ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
          else
            x = cadastruct(x,VarStrings{Vcount},[],0);
          end
        end
      end
      varargout{Vcount} = x;
    end
  elseif strcmp(FunString,'keyboard')
    
  elseif strcmp(FunString,'return')
    
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  if strcmp(FunString,'global')
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if isnumeric(x)
        if all(size(x)) == 0
          % Check to see if need to make a cadastruct
          xid = ADIGATOR.VARINFO.COUNT;
          xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
          asgnlocs = find(ADIGATOR.VARINFO.NAMELOCS(:,1)==xnameloc);
          nextasgn = asgnlocs(asgnlocs>xid);
          if ~ADIGATOR.VARINFO.NAMELOCS(nextasgn,3)
            x = cadastruct([],VarStrings{Vcount},xid,1);
          else
            x = adigatorMakeNumeric(x);
          end
        else
          x = adigatorMakeNumeric(x);
        end
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      else
        x = structParse(x,VarStrings{Vcount},0);
        if isstruct(x)
          if ADIGATOR.VARINFO.COUNT == PreOpCount
            % This guy gets his own variable count
            % Check to see if need to make a cadastruct array
            xid = ADIGATOR.VARINFO.COUNT;
            xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
            asgnlocs = find(ADIGATOR.VARINFO.NAMELOCS(:,1)==xnameloc);
            nextasgn = asgnlocs(asgnlocs>xid);
            if ~isempty(nextasgn) && ~ADIGATOR.VARINFO.NAMELOCS(nextasgn(1),3)
              x = cadastruct([],VarStrings{Vcount},xid,1);
              x = cadaOverMap(x);
            else
              x = cadastruct([],VarStrings{Vcount},xid,0);
            end
            ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
          else
            x = cadastruct(x,VarStrings{Vcount},[],0);
          end
        end
      end
      VarStrings{Vcount} = [VarStrings{Vcount},' '];
      varargout{Vcount} = x;
    end
    fprintf(fid,[indent,'global ',cell2mat(VarStrings),'\n']);
  elseif NUMvars == 1
    % Single Variable Assignment - Either Numeric or Structure/Cell
    x = Variables{1};
    if isnumeric(x)
      % Need to give this its own Operation count and turn it into a
      % Symbolic.
      if all(size(x)) == 0
        % Check to see if need to make a cadastruct
        xid = ADIGATOR.VARINFO.COUNT;
        xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
        asgnlocs = find(ADIGATOR.VARINFO.NAMELOCS(:,1)==xnameloc);
        nextasgn = asgnlocs(asgnlocs>xid);
        if ~ADIGATOR.VARINFO.NAMELOCS(nextasgn,3)
          x = cadastruct([],VarStrings{Vcount},xid,1);
        else
          x = adigatorMakeNumeric(x);
        end
      else
        x = adigatorMakeNumeric(x);
      end
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      % Can just print out whatever is to the right of the =
      if ~ADIGATOR.EMPTYFLAG
        EqLoc = strfind(FunString,'=');
        if isa(x,'cada')
          xfunstr = x.func.name;
          fprintf(fid,[indent,xfunstr,' = ',FunString(EqLoc+1:end),'\n']);
        else
          fprintf(fid,[indent,VarStrings{Vcount},' = [];\n']);
        end
      end
      x = cadaOverMap(x);
    elseif (iscell(x) || isstruct(x))
      % Cell or Structure
      if ~ADIGATOR.EMPTYFLAG
        FunString = FindComments(FunString);
        fprintf(fid,[indent,FunString,'\n']);
      end
      x = structParse(x,VarStrings{Vcount},0);
      if isstruct(x)
        if ADIGATOR.VARINFO.COUNT == PreOpCount
          % This guy gets his own variable count
          % Check to see if need to make a cadastruct array
          xid = ADIGATOR.VARINFO.COUNT;
          xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
          asgnlocs = find(ADIGATOR.VARINFO.NAMELOCS(:,1)==xnameloc);
          nextasgn = asgnlocs(asgnlocs>xid);
          if ~isempty(nextasgn) && ~ADIGATOR.VARINFO.NAMELOCS(nextasgn(1),3)
            x = cadastruct([],VarStrings{Vcount},xid,1);
            x = cadaOverMap(x);
          else
            x = cadastruct([],VarStrings{Vcount},xid,0);
          end
          ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
        else
          x = cadastruct(x,VarStrings{Vcount},[],0);
        end
      end
    elseif ~ADIGATOR.EMPTYFLAG
      FunString = FindComments(FunString);
      fprintf(fid,[indent,FunString,'\n']);
    end
    varargout{Vcount} = x;
  elseif NUMvars > 1
    % Need to build what is to the left of the equals sign before we print
    % out anything. After we print it out, then we can call the OverMap In
    % this case, we should never have cells/structures because we didnt
    % allow in the trace run.
    LHSstrings = cell(1,NUMvars);
    for Vcount = 1:NUMvars
      % x is numeric.
      x = adigatorMakeNumeric(Variables{Vcount});
      if ~ADIGATOR.EMPTYFLAG
        LHSstrings{Vcount} = [x.func.name,','];
      end
      varargout{Vcount}  = x;
    end
    % Print out our LHS and the users RHS
    LHSstrings = cell2mat(LHSstrings);
    if ~ADIGATOR.EMPTYFLAG
      EQloc = strfind(FunString,'=');
      fprintf(fid,[indent,LHSstrings(1:end-1),' = ',...
        FunString(EQloc+1:end),'\n']);
    end
    % Make another run through the variables and check the overmapping.
    for Vcount = 1:NUMvar
      varargout{Vcount} = cadaOverMap(x);
    end
  elseif strcmp(FunString(1),'%')
    % User Comment.
    FunString = FindComments(FunString);
  elseif strcmp(FunString,'keyboard')
    % KEYBOARD
  elseif strcmp(FunString,'return')
    % RETURN STATEMENT
  elseif strfind(FunString,'global')
    
  end

  if ADIGATOR.DERNUMBER == 1 && ADIGATOR.OPTIONS.COMMENTS
    fprintf(fid,[indent,'%%User Line: ',regexprep(FunString,'\\','\\\\'),'\n']);
  elseif length(FunString) > 1 && strcmp(FunString(1),'%')
    fprintf(fid,[indent,regexprep(FunString,'\\','\\\\'),'\n']);
  elseif ADIGATOR.OPTIONS.COMMENTS
    fprintf(fid,[indent,'%% Deriv %1.0d Line: ',...
      regexprep(FunString,'\\','\\\\'),'\n'],ADIGATOR.DERNUMBER-1);
  end
end
ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
ADIGATOR.SUBSINDEXFLAG = 0;
return
end

function FunStri = FindComments(FunStri)

DoinkLocs = strfind(FunStri,'%');
if ~isempty(DoinkLocs)
  NUMdoinks = length(DoinkLocs);
  FunStri(end+NUMdoinks) = ' ';
  for Dcount = 1:length(DoinkLocs)
    Dloc = DoinkLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'%',FunStri(Dloc:end-1)];
  end
end
return
end

function cadaVars = cadaGetGlobalVars(cadaGlobalStrs)

cadaVars = cell(size(cadaGlobalStrs));
for cadaVcount = 1:length(cadaGlobalStrs)
  eval(['global ',cadaGlobalStrs{cadaVcount},';']);
  cadaVars{cadaVcount} = eval(cadaGlobalStrs{cadaVcount});
end

end

function x = structParse(x,xStr,structflag)
global ADIGATOR
if isa(x,'cadastruct') && ~cadaIsArray(x);
  x = cadaGetStruct(x);
end

if isa(x,'cada')
  ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
  if ~structflag
    %xold = x;
    x    = adigatorStructAnalyzer(x,xStr,0);
  end
elseif isnumeric(x)
  if ~structflag
    x = adigatorMakeNumeric(x);
    x = adigatorStructAnalyzer(x,xStr,0);
  elseif ~isempty(x)
    x = adigatorMakeNumeric(x);
  end
elseif isa(x,'cadastruct')
  if structflag
    xID = cadaGetStructID(x);
    ADIGATOR.VARINFO.LASTOCC(xID,1) = ADIGATOR.VARINFO.COUNT;
    x = cadaGetStruct(x);
  else
    x = adigatorStructAnalyzer(x,xStr,0);
  end
elseif isstruct(x)
  x = orderfields(x);
  fnames = fieldnames(x);
  if numel(x) > 1
    for I = 1:numel(x)
      for J = 1:length(fnames)
        F = fnames{J};
        x(I).(F) = structParse(x(I).(F),[xStr,'.',F],1);
      end
    end
    if ~structflag
      x = cadastruct(x,xStr,[],1);
      x = adigatorStructAnalyzer(x,xStr,0);
    end
  elseif ~isempty(x)
    for I = 1:length(fnames)
      F = fnames{I};
      x.(F) = structParse(x.(F),[xStr,'.',F],structflag);
    end
    if ~structflag && ADIGATOR.RUNFLAG > 0
      % Check to see if we need to convert this to a cadastruct array
      lastid = ADIGATOR.VARINFO.COUNT-1;
      lastname = ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(lastid,1)};
      if length(xStr) == length(lastname) && strcmp(lastname,xStr)
        % Convert
        x = cadastruct(x,xStr,lastid,1);
        x = cadaOverMap(x);
      end
    end
  end
elseif iscell(x)
  for I = 1:numel(x)
    x{I} = structParse(x{I},xStr,1);
  end
  if ~structflag
    x = cadastruct(x,xStr,[],1);
    x = adigatorStructAnalyzer(x,xStr,0);
  end
end

end