function varargout = adigatorVarAnalyzer(FunString,varargin)
% CADA overloaded version of adigatorVarAnalyzer
%
% This module is the Overloaded Version of adigatorVarAnalyzer. This is
% called from the intermediate functions after a line of user code has
% been evaluated in order to analyze the outputs. As this is the overloaded
% version, it will only be called if at least one input is overloaded.
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
% changed
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
PreOpCount = ADIGATOR.PREOPCOUNT;

% Parse the inputs;
NUMvars = nargout;
Variables = cell(NUMvars,1);
VarStrings = cell(NUMvars,1);
SubsFlags = zeros(NUMvars,1);
for Vcount = 1:NUMvars
  Variables{Vcount}  = varargin{1+(Vcount-1)*3};
  VarStrings{Vcount} = varargin{2+(Vcount-1)*3};
  SubsFlags(Vcount)  = varargin{3+(Vcount-1)*3};
end
varargout = cell(NUMvars,1);


if ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  for Vcount = 1:NUMvars
    x = Variables{Vcount};
    if PreOpCount == ADIGATOR.VARINFO.COUNT
      % No Overloaded Operation was performed - Direct Assignment
      if isa(x,'cada')
        % Need to give this guy his own operation count
        ADIGATOR.VARINFO.LASTOCC([ADIGATOR.VARINFO.COUNT,...
          x.id],1) = ADIGATOR.VARINFO.COUNT;
        % Set his VarName
        adigatorAssignImpVarNames(ADIGATOR.VARINFO.COUNT,VarStrings{Vcount},SubsFlags(Vcount));
        if isinf(ADIGATOR.VARINFO.NAMELOCS(x.id,3))
          ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = ADIGATOR.VARINFO.NAMELOCS(x.id,3);
        end
        x.id = ADIGATOR.VARINFO.COUNT;
        ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      else
        error(['??? Unable to Process Statement: ',FunString]);
      end
    elseif isa(x,'cada')
      % Just set the Variable Names up and Intermediate Variables.
      xID = x.id;
      adigatorAssignImpVarNames(xID,VarStrings{Vcount},SubsFlags(Vcount));
      % --Set Intermediate Variable Identifiers--
      if xID > PreOpCount && Vcount == 1
        ADIGATOR.VARINFO.NAMELOCS(PreOpCount:xID-1,2)=1:xID-PreOpCount;
      end
      if any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(xID,1)==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
        ADIGATOR.VARINFO.NAMELOCS(xID,3) = -Inf;
      end
    else
      error(['??? Unable to Process Statement: ',FunString]);
    end
    varargout{Vcount} = x;
    if ADIGATOR.DERNUMBER > 1 && ADIGATOR.FORINFO.FLAG && Vcount == 1
      % If in a FOR loop and this doesnt need derivatives printed, make
      % sure that any intermediate operations prior to this one don't get
      % derivatives printed either.
      ADIGATOR.PRINT.FLAG = 1;
      [~,DPFLAG] = cadafuncname(x.id);
      ADIGATOR.PRINT.FLAG = 0;
      if ~DPFLAG
        prevvars = find(ADIGATOR.VARINFO.NAMELOCS(1:end-1,1),1,'last');
        cadaCancelDerivs(x.id,prevvars);
      end
    end
  end
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  for Vcount = 1:NUMvars
    x = Variables{Vcount};
    if PreOpCount == ADIGATOR.VARINFO.COUNT
      % No Overloaded Operation was performed - is a Direct Assignment.
      OverID = x.id;
      x.id   = ADIGATOR.VARINFO.COUNT;
      % OverMapping
      x = cadaOverMap(x,OverID);
      if ~isa(x,'cada')
        x = adigatorMakeNumeric(x);
      end
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
    else
      % Overloaded Operation Performed - Output is Overloaded.
      x = cadaOverMap(x);
    end
    if any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(x.id,1)==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
      derflag = cadaCheckForDerivs(x);
      if derflag && ADIGATOR.DERNUMBER == 1
        error(['variable ''',ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(x.id,1)},...
          ''' is either an auxiliary or global variable which was ',...
          're-assigned to have derivative information - this is not allowed.']);
      elseif derflag
        ADIGATOR.VARINFO.NAMELOCS(x.id,3) = 1;
      end
    end
    varargout{Vcount} = x;
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  if PreOpCount == ADIGATOR.VARINFO.COUNT
    % Direct Assignment - Need to Print this out and to give output proper
    % names.
    x = Variables{1};
    xIDold = x.id;
    xID = ADIGATOR.VARINFO.COUNT;
    [funcname,DPflag] = cadafuncname(xID);
    x.id = xID;
    if ~ADIGATOR.EMPTYFLAG
%       if ADIGATOR.DERNUMBER==1 && isinf(ADIGATOR.VARINFO.NAMELOCS(xIDold,3)) ...
%           && ~isempty(regexp(funcname,'\<cadainput\d*_\d*.f\>','once')) 
%         % Auxiliary data - do not make it a structure
%         funcname = funcname(1:end-2);
%       end
      for Vcount = 1:ADIGATOR.NVAROFDIFF
        if ~isempty(x.deriv(Vcount).nzlocs) 
          if any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(xID,1)==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
            if ADIGATOR.DERNUMBER == 1
              error(['variable ''',ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(x.id,1)},...
                ''' is either an auxiliary or global variable which was ',...
                're-assigned to have derivative information - this is not allowed.']);
            else
              ADIGATOR.VARINFO.NAMELOCS(x.id,3) = 1;
            end
          end
          derivname = cadadername(funcname,Vcount);
          if DPflag
            fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,derivname,' = '...
              x.deriv(Vcount).name,'; ']);
          end
          x.deriv(Vcount).name = derivname;
        end
      end
      fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,funcname,' = ',...
        x.func.name,';\n']);
    end
    x.func.name  = funcname;
    varargout{1} = cadaOverMap(x);
    
    ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  else
    % Need to check and make sure everything is OverMapped properly.
    for Vcount = 1:NUMvars
      x = Variables{Vcount};
      if any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(x.id,1)==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
        derflag = cadaCheckForDerivs(x);
        if derflag 
          if ADIGATOR.DERNUMBER == 1
            error(['variable ''',ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(x.id,1)},...
              ''' is either an auxiliary or global variable which was ',...
              're-assigned to have derivative information - this is not allowed.']);
          else
            ADIGATOR.VARINFO.NAMELOCS(x.id,3) = 1;
          end
        end
      end
      varargout{Vcount} = cadaOverMap(x);
    end
  end
  SlashLocs = strfind(FunString,'\');
  if ~isempty(SlashLocs)
    for Dcount = 1:length(SlashLocs)
      Dloc = SlashLocs(Dcount)+Dcount;
      FunString = [FunString(1:Dloc-1),'\',FunString(Dloc:end-1)];
    end
  end
  if ADIGATOR.OPTIONS.COMMENTS && ~isempty(FunString)
    if ADIGATOR.DERNUMBER == 1 
      fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,'%%User Line: ',...
        FunString,'\n']);
    else
      fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,'%% Deriv %1.0d Line: ',...
        FunString,'\n'],ADIGATOR.DERNUMBER-1);
    end
  end
end
ADIGATOR.PREOPCOUNT = ADIGATOR.VARINFO.COUNT;
ADIGATOR.SUBSINDEXFLAG = 0;
end