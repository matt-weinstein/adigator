function [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs)
% function [FunctionInfo, Outputs] = adigatorFunctionEnd(FunID,FunctionInfo,Outputs)
% This transformation routine is called after the evaluation of an
% intermediate function. It may also be called from within
% adigatorFunctionInitialize if it is determined that we do not need to
% evaluate a particular function.
%
% Inputs:
%   FunID - integer identifying the function
%   FunctionInfo - data structure containing all information on the program
%                  being differentiated
%   Outputs - cell array containing outputs of the function
%
% Outputs:
%   FunctionInfo - data structure containing all information on the program
%                  being differentiated
%   Outputs - cell array containing the (possibly modified) outputs of the 
%             function
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.OPTIONS.PREALLOCATE
  FunctionInfo(FunID).STRUCASGN = ADIGATOR.STRUCASGN;
  
  ADIGATOR.FILE.PARENTID(end) = [];
  if ~isempty(ADIGATOR.FILE.PARENTID)
    ParentID = ADIGATOR.FILE.PARENTID(end);
    ADIGATOR.STRUCASGN = FunctionInfo(ParentID).STRUCASGN;
  end
  return
end
NUMoutputs = length(Outputs);

% ----------------------------------------------------------------------- %
%                             Parse Outputs                               %
% ----------------------------------------------------------------------- %

NumOutVars = 0;
OutVars    = cell(NUMoutputs,1);
OutVarStrs = cell(NUMoutputs,1);
if FunID > 1
  OutOffset = FunctionInfo(ADIGATOR.FILE.PARENTID(end-1)).VARINFO.COUNT-1;
else
  OutOffset = 0;
end

for Ocount = 1:NUMoutputs
  % We are looking for:
  %   1. Overloaded Objects - do nothing with these
  %   2. Numeric Objects - turn these into overloaded objects (according
  %     to what option the user has set, default is that they hold the
  %     same sparsity pattern.
  %   3. Structures/Cells - for these we need to scan them until we get
  %     to either an Overloaded Object or a Numeric Object.
  CurVar    = Outputs{Ocount};
  CurVarStr = FunctionInfo(FunID).Output.Names{Ocount};
  [Outputs{Ocount},OutVars,OutVarStrs,NumOutVars] =...
    FindOutputVars(CurVar,CurVarStr,OutVars,OutVarStrs,NumOutVars,0,OutOffset);
end
OutVars    = OutVars(1:NumOutVars);
OutVarStrs = OutVarStrs(1:NumOutVars);

% OutVars are now all of the overloaded variables which the function starts
% with. OutVarStrs are their string names (with structure and cell reference
% extensions).

if ~ADIGATOR.RUNFLAG
  %% ~~~~~~~~~~~~~~~~~~~~~~~~~ EMPTY RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  FunctionInfo(FunID).Output.StrucNames = OutVarStrs;
  FunctionInfo(FunID).Output.EmptyStruc = OutVars;
  FunctionInfo(FunID).Output.EmptyVars  = Outputs;
  if ~FunctionInfo(FunID).Iteration.CallCount
    % Actually ran this function
    FunctionInfo(FunID).Iteration.CallCount = 1;
    if FunctionInfo(FunID).FunAsLoopFlag
      adigatorForIterEnd(1,0);
    end
    FunctionInfo(FunID).BREAKLOCS = ADIGATOR.BREAKLOCS;
    FunctionInfo(FunID).CONTLOCS  = ADIGATOR.CONTLOCS;
    FunctionInfo(FunID).ERRORLOCS = ADIGATOR.ERRORLOCS;
    CallCount = 1;
  else
    CallCount = 0;
  end
  FunctionInfo = StoreGlobalVars(FunID,FunctionInfo,CallCount);
  if CallCount && FunctionInfo(FunID).DERNUMBER > 1 && ...
      strcmp(FunctionInfo(FunID).PreviousDerivData.Derivative...
      (FunctionInfo(FunID).DERNUMBER-1).FunType,'Main')
    % Take care of Auxiliary Inputs
    PrevAux = FunctionInfo(FunID).PreviousDerivData.AuxInputs;
    if FunID > 1
      PrevAuxFlags = FunctionInfo(FunID).PreviousDerivData.AuxInputFlags;
    else
      PrevAuxFlags = zeros(size(PrevAux));
    end
    InputStrucNames = FunctionInfo(FunID).Input.StrucNames;
    InputCounts = 1:length(InputStrucNames);
    for PAcount = 1:length(PrevAux)
      InputLoc = InputCounts(strcmp(PrevAux{PAcount},InputStrucNames));
      if (FunID > 1 || cadaCheckForDerivs(FunctionInfo(FunID).Input.StrucVars{InputLoc})) &&...
          ~PrevAuxFlags(PAcount)
        FunctionInfo(FunID).VARINFO.NAMES{InputLoc} = ...
          [FunctionInfo(FunID).VARINFO.NAMES{InputLoc},'.f'];
      end
    end
  end
  if FunID > 1
    DPflag = FunctionInfo(FunID).DERNUMBER > 1 && ADIGATOR.DERNUMBER == 1;
    % ------- Assign Naming Scheme for Output in Calling Func ----------- %
    OutputStrs = FunctionInfo(FunID).Output.Names.';
    NewNames   = cell(1,NUMoutputs);
    NameChecks = cell(1,NUMoutputs);
    for Ocount = 1:NUMoutputs
      NewNames{Ocount} = sprintf('cadaoutput%1.0f_%1.0f',FunID,Ocount);
      NameChecks{Ocount} = ['\<',OutputStrs{Ocount},'\W'];
    end
    for Ocount = 1:NumOutVars
      OutVarStr = OutVarStrs{Ocount};
      Matches = regexp(OutVarStr,NameChecks,'once','end');
      for Mcount = 1:length(Matches);
        if isempty(Matches{Mcount}); Matches{Mcount} = 0; end;
      end
      Matches = cell2mat(Matches);
      if any(Matches)
        OutStr  = [NewNames{logical(Matches)},OutVarStr(nonzeros(Matches):end)];
      else
        Matches = strcmp(OutVarStr,OutputStrs);
        OutStr  = NewNames{Matches};
      end
      if isa(OutVars{Ocount},'cada')
        varID = OutVars{Ocount}.id;
      else
        varID = cadaGetStructID(OutVars{Ocount});
      end
      
      adigatorAssignImpVarNames(varID,OutStr,0);
      if DPflag
        ADIGATOR.VARINFO.NAMELOCS(varID,2) = -FunID;
      end
    end
  end
elseif ADIGATOR.OPTIONS.UNROLL && FunID > 1
  %% ~~~~~~~~~~~~~~~~~~~ UNROLLING SUB-FUNCTION ~~~~~~~~~~~~~~~~~~~~~~~~ %%
  if ADIGATOR.RUNFLAG == 2
    CallCount = FunctionInfo(FunID).Iteration.CallCount;
    IterCount = FunctionInfo(FunID).Iteration.IterID(CallCount);
    FunctionInfo(FunID).Output.Vars(:,IterCount)      = Outputs;
    FunctionInfo(FunID).Output.StrucVars(:,IterCount) = OutVars;
    SaveDataFile(FunID,FunctionInfo);
    fprintf(ADIGATOR.PRINT.FID,'end\n');
  end
  
  FunctionInfo = StoreGlobalVars(FunID,FunctionInfo,1);
  ADIGATOR.PRINT.FID = FunctionInfo(FunID).Iteration.CallerID(1);
  ADIGATOR.RUNFLAG   = FunctionInfo(FunID).Iteration.CallerID(2);
  if ADIGATOR.RUNFLAG == 2
    for Ocount = 1:NumOutVars
      x = OutVars{Ocount};
      if isa(x,'cada')
        xID = x.id;
        funcstr = cadafuncname(xID);
        x.func.name = funcstr;
        for Vcount = 1:ADIGATOR.NVAROFDIFF
          if ~isempty(x.deriv(Vcount).nzlocs)
            x.deriv(Vcount).name = cadadername(funcstr,Vcount,xID);
          end
        end
      else
        [xv,~,xid,xflag] = cadastructDecomp(x);
        xname = ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(xid,1)};
        x = cadastruct(xv,xname,xid,xflag);
      end
      OutVars{Ocount} = x;
    end
    Outputs = cadaAssignOutputs(FunctionInfo(FunID).Output.Names,OutVarStrs,OutVars);
  end
elseif ADIGATOR.RUNFLAG == 1
  %% ~~~~~~~~~~~~~~~~~~~~~~~~~ OVERMAP RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  if ~ADIGATOR.EMPTYFLAG
    CallCount = FunctionInfo(FunID).Iteration.CallCount;
    IterCount = FunctionInfo(FunID).Iteration.IterID(CallCount);
    FunctionInfo(FunID).Output.Vars(:,IterCount)      = Outputs;
    FunctionInfo(FunID).Output.StrucVars(:,IterCount) = OutVars;
    SameFlag = any(FunctionInfo(FunID).Iteration.IterID(1:end-1) == FunctionInfo(FunID).Iteration.IterID(end));
    if FunctionInfo(FunID).FunAsLoopFlag && ~SameFlag
      adigatorForIterEnd(1,IterCount);
    end
    FunctionInfo = StoreGlobalVars(FunID,FunctionInfo,IterCount);
  else
    FunctionInfo = StoreGlobalVars(FunID,FunctionInfo,0);
  end
  
  if ADIGATOR.FORINFO.FLAG
    for Ocount = 1:NumOutVars
      x = OutVars{Ocount};
      if isa(x,'cada')
        xID = x.id;
        x = cadaOverMap(x);
        x.id = xID;
      else
        xID = cadaGetStructID(x);
        x = cadaOvermap(x);
        [xv,xname,~,xflag] = cadastructDecomp(x);
        x = cadastruct(xv,xname,xID,xflag);
      end
      OutVars{Ocount}  = x;
    end
  end
else
  %% ~~~~~~~~~~~~~~~~~~~~~~~~ PRINTING RUN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %%
  if ~FunctionInfo(FunID).Iteration.CallCount && ~ADIGATOR.EMPTYFLAG
    % ------------------- Printing This Function ------------------------ %
    if FunID == 1
      for Ocount = 1:NumOutVars
        x = OutVars{Ocount};
        if isa(x,'cada')
          adigatorPrintOutputIndices(x);
        elseif isa(x,'cadastruct')
          adigatorPrintOutputIndices(x);
          x = cadaGetStruct(x);
        end
        OutVars{Ocount} = x;
      end
      if FunctionInfo(FunID).DERNUMBER > 1
        
      end
    end
    SaveDataFile(FunID,FunctionInfo);
    fprintf(ADIGATOR.PRINT.FID,'end\n');
    Outputs = cadaAssignOutputs(FunctionInfo(FunID).Output.Names,OutVarStrs,OutVars);
  elseif ADIGATOR.EMPTYFLAG
    FunctionInfo = StoreGlobalVars(FunID,FunctionInfo,0);
  else
    % ------------------- Printing Parent Function ---------------------- %
    
    FunctionInfo = StoreGlobalVars(FunID,FunctionInfo,0);
    DPflag = FunctionInfo(FunID).DERNUMBER > 1 && ADIGATOR.DERNUMBER == 1;
    
    if ADIGATOR.FORINFO.FLAG
      % Parent Function is in a FOR loop - Get outputs from VARSTORAGE
      ParentID = ADIGATOR.FILE.FUNID;
      for Ocount = 1:NumOutVars
        x   = OutVars{Ocount};
        if isa(x,'cada')
          xID = x.id;
          OverLoc            = ADIGATOR.VARINFO.OVERMAP.FOR(xID,1);
          OutVars{Ocount}    = FunctionInfo(ParentID).VARSTORAGE.OVERMAP{OverLoc};
          OutVars{Ocount}.id = xID;
        else
          xID = cadagetStructID(x);
          OverLoc            = ADIGATOR.VARINFO.OVERMAP.FOR(xID,1);
          x = FunctionInfo(ParentID).VARSTORAGE.OVERMAP{OverLoc};
          [xv,xname,~,xflag] = cadastructDecomp(x);
          OutVars{Ocount}    = cadastruct(xv,xname,xID,xflag);
        end
      end
    end
    
    % OutVars are now what should be returned to the calling function - we
    % now need to remap the outputs of the child function to what should
    % return to the parent function. We get the Child Function Overmap from
    % VARSTORAGE.
    % The .id field in OutVars is proper for the Parent function, and
    % cadaPrintRemap will do the naming using this ID
    
    VARINFO = FunctionInfo(FunID).VARINFO;
    
    NumIDs  = VARINFO.COUNT-1;
    AllVarIDs   = 1:NumIDs;
    OutputIDs   = AllVarIDs(VARINFO.LASTOCC == NumIDs+1);
    OutputNames = VARINFO.NAMES(VARINFO.NAMELOCS(OutputIDs,1));
    
    
    for Ocount = 1:NumOutVars
      x     = OutVars{Ocount};
      if isa(x,'cada')
        varID = x.id;
      else
        varID = cadaGetStructID(x);
      end
      
      if FunctionInfo(FunID).FunAsLoopFlag
        OutputID  = OutputIDs(strcmp(OutVarStrs{Ocount},OutputNames));
        OverLoc   = VARINFO.OVERMAP.FOR(OutputID,1);
        OverVar   = FunctionInfo(FunID).VARSTORAGE.OVERMAP{OverLoc};
        x         = cadaPrintReMap(OverVar,x,varID);
      elseif isa(x,'cada')
        funcstr = cadafuncname(varID);
        x.func.name = funcstr;
        for Vcount = 1:ADIGATOR.NVAROFDIFF
          if ~isempty(x.deriv(Vcount).nzlocs)
            x.deriv(Vcount).name = cadadername(funcstr,Vcount,varID);
          end
        end
      else
        [x,~,xid,xflag] = cadastructDecomp(x);
        xname = ADIGATOR.VARINFO.NAMES{ADIGATOR.VARINFO.NAMELOCS(xid,1)};
        x = cadastruct(x,xname,xid,xflag);
      end
      OutVars{Ocount} = x;
    end
    
    % Need to now assign OutVars to Outputs.
    Outputs = cadaAssignOutputs(FunctionInfo(FunID).Output.Names,OutVarStrs,OutVars);
  end
end

ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+NumOutVars;
ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
end

function [CurVar,Vars,VarStrs,NumOutVars] =...
  FindOutputVars(CurVar,CurVarStr,Vars,VarStrs,NumOutVars,NameFlag,OutOffset)
if isa(CurVar,'cadastruct') && ~cadaIsArray(CurVar)
  CurVar = cadaGetStruct(CurVar);
end
global ADIGATOR
if isa(CurVar,'cada')
  NumOutVars = NumOutVars+1;
  if length(Vars) < NumOutVars
    Vars{NumOutVars*2,1}    = [];
    VarStrs{NumOutVars*2,1} = [];
  end
  
  if NameFlag
    funcstr = cadafuncname(CurVar.id);
    CurVar.func.name = funcstr;
    for Vcount = 1:ADIGATOR.NVAROFDIFF
      if ~isempty(CurVar.deriv(Vcount).nzlocs)
        CurVar.deriv(Vcount).name = cadadername(funcstr,Vcount,CurVar.id);
      end
    end
  else
    ADIGATOR.VARINFO.LASTOCC(CurVar.id,1) = size(ADIGATOR.VARINFO.NAMELOCS,1)+1;
    if OutOffset
      CurVar.id = NumOutVars+OutOffset;
    end
  end
  Vars{NumOutVars}    = CurVar;
  VarStrs{NumOutVars} = CurVarStr;
elseif isa(CurVar,'cadastruct')
  % Structure Array.
  NumOutVars = NumOutVars+1;
  if length(Vars) < NumOutVars
    Vars{NumOutVars*2,1}    = [];
    VarStrs{NumOutVars*2,1} = [];
  end
  
  if ~NameFlag
    [xval, xname, xid] = cadastructDecomp(CurVar);
    ADIGATOR.VARINFO.LASTOCC(xid,1) = size(ADIGATOR.VARINFO.NAMELOCS,1)+1;
    if OutOffset
      CurVar = cadastruct(xval,xname,NumOutVars+OutOffset,1);
    end
  end
  Vars{NumOutVars}    = CurVar;
  VarStrs{NumOutVars} = CurVarStr;
elseif isstruct(CurVar) && ~isempty(CurVar)
  % Scalar Structure
  CurVar = orderfields(CurVar);
  StructNames = fieldnames(CurVar);
  for Scount = 1:length(StructNames)
    NewVar    = CurVar.(StructNames{Scount});
    NewVarStr = [CurVarStr,'.',StructNames{Scount}];
    [CurVar.(StructNames{Scount}),Vars,VarStrs,NumOutVars] =...
      FindOutputVars(NewVar,NewVarStr,Vars,VarStrs,NumOutVars,NameFlag,OutOffset);
  end
else
  % Do nothing.
end
end

function FunctionInfo = StoreGlobalVars(FunID, FunctionInfo,CallCount)
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE ADIGATORDATA

if CallCount
  FunctionInfo(FunID).VARINFO    = ADIGATOR.VARINFO;
  FunctionInfo(FunID).FORDATA    = ADIGATORFORDATA;
  FunctionInfo(FunID).IFDATA     = ADIGATOR.IFDATA;
  FunctionInfo(FunID).FORINFO    = ADIGATOR.FORINFO;
  FunctionInfo(FunID).VARSTORAGE = ADIGATORVARIABLESTORAGE;
  FunctionInfo(FunID).STRUCASGN = ADIGATOR.STRUCASGN;
end

if FunID > 1
  ADIGATOR.FILE.PARENTID(end) = [];
  ParentID = ADIGATOR.FILE.PARENTID(end);
  ADIGATOR.FILE.FUNID = ParentID;
  
  ADIGATOR.STRUCASGN = FunctionInfo(ParentID).STRUCASGN;
  ADIGATORVARIABLESTORAGE = FunctionInfo(ParentID).VARSTORAGE;
  ADIGATORFORDATA         = FunctionInfo(ParentID).FORDATA;
  ADIGATOR.IFDATA         = FunctionInfo(ParentID).IFDATA;
  ADIGATOR.FORINFO        = FunctionInfo(ParentID).FORINFO;
  ADIGATOR.VARINFO        = FunctionInfo(ParentID).VARINFO;
  ADIGATOR.DERNUMBER      = FunctionInfo(ParentID).DERNUMBER;
  if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.RUNFLAG
    ADIGATORDATA = FunctionInfo(ParentID).DATA;
  end
end
end

function SaveDataFile(FunID,FunctionInfo)
global ADIGATOR ADIGATORDATA
ADiGatorMatFileName = ADIGATOR.PRINT.FILENAME;
if ADIGATOR.DERNUMBER > 1
  PrevDerivStruc = FunctionInfo(FunID).PreviousDerivData;
  ADiGatorstruc.Derivative = PrevDerivStruc.Derivative;
  for Dcount = 1:ADIGATOR.DERNUMBER-1
    %IndStr = sprintf('Gator%1.0dIndices',Dcount);
    %ADiGatorstruc.(IndStr) = PrevDerivStruc.(IndStr);
    DatStr = sprintf('Gator%1.0dData',Dcount);
    ADiGatorstruc.(DatStr) = PrevDerivStruc.(DatStr);
  end
end
if FunID == 1
  ADiGatorFunName = ADIGATOR.PRINT.FILENAME;
  ADiGatorstruc.Derivative(ADIGATOR.DERNUMBER).FunType = 'Main';
else
  if FunctionInfo(FunID).DERNUMBER > 1
    ADiGatorFunName = FunctionInfo(FunID).File.Name;
    if length(ADiGatorFunName) > 8 && strcmp(ADiGatorFunName(1:7),'ADiGator')
      ADiGatorFunName = ['ADiGator_',ADiGatorFunName(8:end)];
    end
  elseif ADIGATOR.OPTIONS.UNROLL && FunID > 1
    IterCount = FunctionInfo(FunID).Iteration.IterCount;
    ADiGatorFunName = sprintf(['ADiGator_',FunctionInfo(FunID).File.Name,'%1.0f'],IterCount);
  else
    ADiGatorFunName = ['ADiGator_',FunctionInfo(FunID).File.Name];
  end
  ADiGatorstruc.Derivative(ADIGATOR.DERNUMBER).FunType = 'Sub';
end
ADiGatorstruc.Derivative(ADIGATOR.DERNUMBER).Variables = ADIGATOR.VAROFDIFF;
if FunID == 1
  % Build some input checks for derivatives
  AuxInputs     = FunctionInfo(FunID).Input.StrucNames;
  CurrentInputs = FunctionInfo(FunID).Input.StrucVars;
  InputChecks   = struct('name',[],'deriv',[]);
  for Icount = 1:length(FunctionInfo(FunID).Input.StrucVars)
    curInput = CurrentInputs{Icount};
    flag = cadaCheckForDerivs(curInput);
    if flag
      AuxInputs{Icount} = [];
    end
  end
  ADiGatorstruc.Derivative(ADIGATOR.DERNUMBER).InputChecks = InputChecks;
  if ADIGATOR.DERNUMBER > 1
    OldAuxInputs = FunctionInfo(1).PreviousDerivData.AuxInputs;
    for Acount = 1:length(AuxInputs)
      if ~any(strcmp(AuxInputs{Acount},OldAuxInputs))
        AuxInputs{Acount} = [];
      end
    end
  end
  ADiGatorstruc.AuxInputs = AuxInputs(~cellfun(@isempty,AuxInputs));
end

%ADiGatorstruc.(sprintf('Gator%1.0dIndices',ADIGATOR.DERNUMBER)) = ADIGATORDATA.INDICES;
ADiGatorstruc.(sprintf('Gator%1.0dData',ADIGATOR.DERNUMBER)) = ADIGATORDATA.DATA; %#ok<STRNU>

ADiGatorCallingDir = cd;
eval([ADiGatorFunName,' = ADiGatorstruc;']);
if ~exist([ADiGatorCallingDir,filesep,ADiGatorMatFileName,'.mat'],'file');
  save([ADiGatorCallingDir,filesep,ADiGatorMatFileName,'.mat'],ADiGatorFunName);
else
  save([ADiGatorCallingDir,filesep,ADiGatorMatFileName,'.mat'],ADiGatorFunName,'-append');
end
rehash
end

function cadaOutVars = cadaAssignOutputs(cadaOutputStrs,cadaOutputStrsFull,cadaOutVarsFull) %#ok<INUSD>
cadaOutVars = cell(size(cadaOutputStrs));
for cadaOcount = 1:length(cadaOutputStrsFull)
  eval([cadaOutputStrsFull{cadaOcount},' = cadaOutVarsFull{cadaOcount};']);
end
for cadaOcount = 1:length(cadaOutVars)
  cadaOutVars{cadaOcount} = eval(cadaOutputStrs{cadaOcount});
end
end