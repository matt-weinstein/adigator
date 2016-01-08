function x = adigatorVarAnalyzer(FunString,x,xStr,subsflag)
%function varargout = adigatorVarAnalyzer(FunString,x,xStr,subsflag) (overloaded
%cadstruct version)
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

% NOTE: If subsflag && ~x.arrayflag --> this is all taken care of within
% overloaded subsasgn.


global ADIGATOR
PreOpCount = ADIGATOR.PREOPCOUNT;
fid     =   ADIGATOR.PRINT.FID;
PFLAG   =   ADIGATOR.PRINT.FLAG;
indent  =   ADIGATOR.PRINT.INDENT;
if ADIGATOR.FILE.CALLFLAG
  % Function was called
  PreOpCount = adigatorFunctionCalled(Variables,VarStrings,SubsFlags);
  if ~PreOpCount
    % adigatorFunctionCalled took care of everything.
    ADIGATOR.PREOPCOUNT = ADIGATOR.VARINFO.COUNT;
    ADIGATOR.FILE.CALLFLAG = 0;
    if ADIGATOR.PRINT.FLAG
      if ADIGATOR.DERNUMBER == 1
        fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,'%%User Line: ',...
          FunString,'\n\n']);
      else
        fprintf(ADIGATOR.PRINT.FID,[ADIGATOR.PRINT.INDENT,'%% Deriv %1.0d Line: ',...
          FunString,'\n'],ADIGATOR.DERNUMBER-1);
      end
    end
    return
  end
  % Only way this wont fire as a return is that operations were done after
  % the function call. Ex. y = myfunc(x)*a;
end

if x.arrayflag
  if ~ADIGATOR.RUNFLAG
    % ------------------------------------------------------------------- %
    %                            Empty Run                                %
    % ------------------------------------------------------------------- %
    if PreOpCount == ADIGATOR.VARINFO.COUNT
      % Need to give this guy his own operation count
      xID = x.id;
      ADIGATOR.VARINFO.LASTOCC([ADIGATOR.VARINFO.COUNT,...
        xID],1) = ADIGATOR.VARINFO.COUNT;
      % Set his VarName
      adigatorAssignImpVarNames(ADIGATOR.VARINFO.COUNT,xStr,subsflag);
      if isinf(ADIGATOR.VARINFO.NAMELOCS(xID,3))
        ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3) = ADIGATOR.VARINFO.NAMELOCS(xID,3);
      end
      if ADIGATOR.VARINFO.NAMELOCS(xID,2) < 0
        ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,2) = ADIGATOR.VARINFO.NAMELOCS(xID,2);
      end
      x.id   = ADIGATOR.VARINFO.COUNT;
      x.name = xStr;
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
    else
      xID = x.id;
      adigatorAssignImpVarNames(xID,xStr,subsflag);
      
      % --Set Intermediate Variable Identifiers--
      if xID > PreOpCount
        ADIGATOR.VARINFO.NAMELOCS(PreOpCount:xID-1,2)=1:xID-PreOpCount;
      end
      if any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(xID,1)==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
        ADIGATOR.VARINFO.NAMELOCS(xID,3) = -Inf;
      end
    end
  elseif ADIGATOR.RUNFLAG == 1
    % ------------------------------------------------------------------- %
    %                          Overmap Run                                %
    % ------------------------------------------------------------------- %
    if PreOpCount == ADIGATOR.VARINFO.COUNT
      % No Overloaded Operation was performed - is a Direct Assignment.
      OverID = x.id;
      x.id   = ADIGATOR.VARINFO.COUNT;
      % OverMapping
      x.name = xStr;
      x = cadaOverMap(x,OverID);
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
    else
      x = cadaOverMap(x);
      x.name = xStr;
    end
  else
    % ------------------------------------------------------------------- %
    %                          Printing Run                               %
    % ------------------------------------------------------------------- %
    if PreOpCount == ADIGATOR.VARINFO.COUNT
      % Need to give this guy his own operation count
      oldID = x.id;
      ADIGATOR.VARINFO.LASTOCC([ADIGATOR.VARINFO.COUNT,...
        oldID],1) = ADIGATOR.VARINFO.COUNT;
      x.id   = ADIGATOR.VARINFO.COUNT;
      oldName = x.name;
      if ~ADIGATOR.EMPTYFLAG
        fprintf(fid,[indent,xStr,' = ',oldName,';\n']);
      end
      x.name = xStr;
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
      if any(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.NAMELOCS(x.id,1)==ADIGATOR.VARINFO.NAMELOCS(:,1),3) == -Inf)
        derflag = cadaCheckForDerivs(x);
        if derflag
          error(['variable ''',ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(x.id,1)},...
            ''' is either an auxillory or global variable which was ',...
            're-assigned to have derivative information - this is not allowed.']);
        end
      end
      xid = ADIGATOR.VARINFO.COUNT-1;
      if ~ADIGATOR.EMPTYFLAG && any(ADIGATOR.VARINFO.NAMELOCS(xid,2) ~= ...
          ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.LASTOCC(:,1)==xid,2))
        adigatorPrintStructAsgn(x,xStr,oldName,xid,oldID);
      end
    end
    
    x = cadaOverMap(x);
  end
elseif ~subsflag
  % scalar structure, direct assignment
  if PreOpCount == ADIGATOR.VARINFO.COUNT
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      if strcmp(x.name,'adigatordummystruct')
        % This resulted from cadastruct subsref - dont want to print
        % FunString
        if ~strcmp(xStr,'adigatordummystruct')
          fprintf(fid,[indent,xStr,' = ',x.name,';\n']);
        else
          FunString = '';
        end
      else
        fprintf(fid,[indent,FunString,'\n']);
      end
    end
    x = adigatorStructAnalyzer(x,xStr,0);
    if isstruct(x)
      x = cadastruct(x,xStr,[],0);
    end
  end
end

if ADIGATOR.RUNFLAG == 2
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
ADIGATOR.FILE.CALLFLAG = 0;
ADIGATOR.SUBSINDEXFLAG = 0;
return
end

