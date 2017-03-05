function [outEvalStr,outEvalVar] = adigatorIfIterEnd(IfCount,BroCount)
% This transformation routine is placed at the end of each IF/ELSEIF/ELSE
% block in the intermediate program.
%
% Inputs: 
%   IfCount  - integer identifying the conditional if/elseif/else set
%   BroCount - integer identifying the branch of the conditional set 
%              (e.g. 1 corresponds to the opening "if" branch)
%
% Outputs:
%   outEvalStr - cell arary of strings to be evaluated on the output in 
%                order to modify the workspace
%   outEvalVar - cell array containing variables to be placed into the 
%                workspace 
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR
if ADIGATOR.OPTIONS.PREALLOCATE
  % Pre-Allocating cells/structures - have to treat everything as if it is
  % a true statement.
  outEvalStr = [];
  outEvalVar = [];
elseif ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  EndCount = ADIGATOR.VARINFO.COUNT-1;
  ADIGATOR.IFDATA(IfCount).BROS(BroCount).END = EndCount;
  
  if size(ADIGATOR.VARINFO.NAMELOCS,1) == EndCount && ~ADIGATOR.FORINFO.FLAG
    % Check for variables which are defined within this branch, and after
    % the last assignment in the branch, are used later within this branch
    StartCount  = ADIGATOR.IFDATA(IfCount).BROS(BroCount).START;
    InnerCounts = StartCount:EndCount;
    InnerNameLocs = nonzeros(unique(ADIGATOR.VARINFO.NAMELOCS(InnerCounts,1)));
    nInnerNameLocs = length(InnerNameLocs);
    INNERDEP = zeros(1,nInnerNameLocs);
    for Icount = 1:nInnerNameLocs
      NameLoc = InnerNameLocs(Icount);
      LastAsgnCount = InnerCounts(find(ADIGATOR.VARINFO.NAMELOCS(InnerCounts,1) == NameLoc,1,'last'));
      if ADIGATOR.VARINFO.LASTOCC(LastAsgnCount) > LastAsgnCount
        % This is a variable which is used after its last assignment in the
        % branch
        INNERDEP(Icount) = LastAsgnCount;
      end
    end
    ADIGATOR.IFDATA(IfCount).BROS(BroCount).INNERDEP = nonzeros(INNERDEP).';
  end
  

  if BroCount > 1 && size(ADIGATOR.VARINFO.NAMELOCS,1) == EndCount
    % Check for Variables which this branch is dependent upon, but are
    % defined in a prior branch.
    IfStart     = ADIGATOR.IFDATA(IfCount).BROS(1).START;
    StartCount  = ADIGATOR.IFDATA(IfCount).BROS(BroCount).START;
    PriorCounts = IfStart:StartCount-1;
    SaveCounts  = PriorCounts(ADIGATOR.VARINFO.LASTOCC(PriorCounts,1) >= StartCount);
    if ~isempty(SaveCounts)
      % This ELSEIF/ELSE section is dependent upon variables changed from
      % within a prior brother IF/ELSEIF statement, will have to store said
      % variables prior to the evaluation of the IF statement, and load
      % them back in prior to the evaluation of this section of code.
      NameLocs   = ADIGATOR.VARINFO.NAMELOCS(SaveCounts,1);
      PrevCounts = 1:IfStart-1;
      DepCounts  = zeros(1,length(NameLocs));
      for Dcount = 1:length(NameLocs)
        deploc = PrevCounts(find(ADIGATOR.VARINFO.NAMELOCS(PrevCounts,1)...
          == NameLocs(Dcount),1,'last'));
        if ~isempty(deploc)
          DepCounts(Dcount) = deploc;
          ADIGATOR.VARINFO.LASTOCC(deploc)  = StartCount;
        else
          error(['Subs-assigning to variable: ''',ADIGATOR.VARINFO.NAMES{NameLocs(Dcount)},...
            ''' within a conditional fragment, without pre-allocating prior to fragment']);
        end
      end
     
      ADIGATOR.VARINFO.LASTOCC(SaveCounts) = StartCount-1;
      ADIGATOR.IFDATA(IfCount).BROS(BroCount).PRIORDEP = DepCounts;
    end
  end
  
  if BroCount == length(ADIGATOR.IFDATA(IfCount).BROS)...
      && ADIGATOR.OPTIONS.UNROLL && ADIGATOR.IFDATA(IfCount).OUTERFLAG...
      && size(ADIGATOR.VARINFO.NAMELOCS,1) == EndCount
    % We are going to run this block twice in a row (once for overmapping
    % and once for printing), need to find variables which occur prior to
    % block, are used within block, and are redefined within block.
    IfStart     = ADIGATOR.IFDATA(IfCount).BROS(1).START;
    PrevCounts  = 1:IfStart-1;
    UsedCounts  = PrevCounts(ADIGATOR.VARINFO.LASTOCC(PrevCounts,1) >= StartCount);
    if ~isempty(UsedCounts)
      UsedNameLocs = ADIGATOR.VARINFO.NAMELOCS(UsedCounts,1);
      IfCounts = StartCount:EndCount;
      SaveCounts = zeros(size(UsedCounts));
      for Scount = 1:length(SaveCounts)
        if any(ADIGATOR.VARINFO.NAMELOCS(IfCounts,1) == UsedNameLocs(Scount))
          SaveCounts(Scount) = UsedCounts(Scount);
        end
      end
      ADIGATOR.IFDATA(IfCount).BROS(1).PRIORDEP = nonzeros(SaveCounts).';
    end
  end
  if BroCount == length(ADIGATOR.IFDATA(IfCount).BROS)
    outEvalStr = [];
    outEvalVar = [];
  end
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  if BroCount == length(ADIGATOR.IFDATA(IfCount).BROS) 
    % Collect any children overmaps and build outEvalVar/outEvalStr such
    % that the overmaps will be assigned to the evaluating workspace if
    % required.
    [outEvalStr outEvalVar] = Returns_InnerUnions_OuterSaves(IfCount,0);
    ADIGATOR.EMPTYFLAG = ADIGATOR.IFDATA(IfCount).EMPTYFLAG;
  else
    outEvalStr = [];
    outEvalVar = [];
  end
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  indent(1:4) = [];
  
  if ADIGATOR.IFDATA(IfCount).BROS(BroCount).RUNFLAG && ...
      ~isempty(ADIGATOR.IFDATA(IfCount).BROS(BroCount).REMAPS)
    % Do any re-mapping that is required
    DoBroRemapping(IfCount,BroCount);
  end
  if BroCount == length(ADIGATOR.IFDATA(IfCount).BROS)
    % Check to see if at least one of the statements ran
    if ~ADIGATOR.FORINFO.FLAG && ADIGATOR.IFDATA(IfCount).PRINTFLAG && ~ADIGATOR.IFDATA(IfCount).ELSEFLAG
      ElseRemapFlag = 1;
      fprintf(fid,[indent,'else\n']);
    else
      ElseRemapFlag = 0;
    end
    
    [outEvalStr outEvalVar] = ...
      Returns_InnerUnions_OuterSaves(IfCount,ElseRemapFlag);

    if ADIGATOR.IFDATA(IfCount).PRINTFLAG
      fprintf(fid,[indent,'end\n']);
    end
    
    ADIGATOR.EMPTYFLAG = ADIGATOR.IFDATA(IfCount).EMPTYFLAG;
  else
    outEvalStr = [];
    outEvalVar = [];
  end
  ADIGATOR.PRINT.INDENT = indent;
end
if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.RUNFLAG &&...
    BroCount == length(ADIGATOR.IFDATA(IfCount).BROS) && ...
    ~isempty(ADIGATOR.IFDATA(IfCount).RESETLOC)
  ADIGATOR.IFINFO.INNERLOC = ADIGATOR.IFDATA(IfCount).RESETLOC;
end
end

function [outEvalStr, outEvalVar] = Returns_InnerUnions_OuterSaves(IfCount,ElseRemapFlag)
global ADIGATOR ADIGATORVARIABLESTORAGE
IfOvermap  = ADIGATOR.IFDATA(IfCount).OVERMAP;
nOvermaps  = length(IfOvermap);
outEvalStr = cell(nOvermaps,1);
outEvalVar = cell(nOvermaps,1);

if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.FORINFO.INNERLOC
  PriorDepStruc = ADIGATOR.IFDATA(IfCount).PRIORDEP;
end

for Ocount = 1:nOvermaps
  OuterOverLoc = IfOvermap(Ocount).Outer(1);
  OuterSaveLoc = IfOvermap(Ocount).Outer(2);
  OuterNameLoc = IfOvermap(Ocount).Return;
  OuterOverObj = ADIGATORVARIABLESTORAGE.OVERMAP{OuterOverLoc};
  PriorSaveLoc  = IfOvermap(Ocount).Prior;
  if ADIGATOR.RUNFLAG == 1 && ~isempty(IfOvermap(Ocount).Inner)
    % Gather any overmaps which belong to this one from children
    % conditional blocks
    for InnerOverLoc = IfOvermap(Ocount).Inner
      InnerOverObj = ADIGATORVARIABLESTORAGE.OVERMAP{InnerOverLoc};
      OuterOverObj = cadaUnionVars(OuterOverObj,InnerOverObj);
    end

    if ~isempty(ADIGATOR.VARINFO.SAVE.FOR)
      if isa(OuterOverObj,'cadastruct')
        outerID = cadaGetStructID(OuterOverObj);
      else
        outerID = OuterOverObj.id;
      end
      ForSaveLoc = ADIGATOR.VARINFO.SAVE.FOR(outerID,2);
      if ForSaveLoc
        ADIGATORVARIABLESTORAGE.SAVE{ForSaveLoc} = OuterOverObj;
      end
    end
    ADIGATORVARIABLESTORAGE.OVERMAP{OuterOverLoc} = OuterOverObj;
  end
  
  if ADIGATOR.RUNFLAG == 2 && OuterSaveLoc && ~ADIGATOR.FORINFO.FLAG
    % Save this overmap so that a parent brother IfIterEnd can do the remap
    ADIGATORVARIABLESTORAGE.SAVE{OuterSaveLoc} = OuterOverObj;
  end
  
  if ~isempty(OuterNameLoc) && (ADIGATOR.RUNFLAG == 1 || ~ADIGATOR.FORINFO.FLAG)
    % Need to return this overmap to the workspace
    VarName = ADIGATOR.VARINFO.NAMES{OuterNameLoc};
    if ADIGATOR.RUNFLAG == 2
      if isa(OuterOverObj,'cadastruct')
        % Change name of cadastruct
        [oVal,~,oID,oAF] = cadastructDecomp(OuterOverObj);
        OuterOverObj = cadastruct(oVal,VarName,oID,oAF);
      else
        % Change name of cada
        outerID = OuterOverObj.id;
        funcname = cadafuncname(outerID);
        OuterOverObj.func.name = funcname;
        for Vcount = 1:ADIGATOR.NVAROFDIFF
          if ~isempty(OuterOverObj.deriv(Vcount).nzlocs)
            OuterOverObj.deriv(Vcount).name = cadadername(funcname,Vcount,outerID);
          end
        end
      end
    end
    outEvalStr{Ocount} = VarName;
    outEvalVar{Ocount} = OuterOverObj;
  end
  
  if ElseRemapFlag && ~isempty(PriorSaveLoc)
    % Must do the Remapping within an imposed ELSE branch
    if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.FORINFO.INNERLOC && ...
        any(PriorDepStruc.Locs == PriorSaveLoc)
      PriorObj = PriorDepStruc.Vars{PriorDepStruc.Locs == PriorSaveLoc};
    else
      PriorObj = ADIGATORVARIABLESTORAGE.SAVE{ADIGATOR.VARINFO.SAVE.IF(PriorSaveLoc,1)};
    end
    if isa(PriorObj,'cadastruct')
      priorID = cadaGetStructID(PriorObj);
    else
      priorID = PriorObj.id;
    end
    cadaPrintReMap(PriorObj,OuterOverObj,priorID);
  end
end
outLogical = ~cellfun(@isempty,outEvalStr);
outEvalStr = outEvalStr(outLogical);
outEvalVar = outEvalVar(outLogical);
for Ocount = 1:length(outEvalStr)
  outEvalStr{Ocount} = sprintf([outEvalStr{Ocount},' = adigatorIfEvalVar{%1.0f};'],Ocount);
end

end



function DoBroRemapping(IfCount,BroCount)
global ADIGATOR ADIGATORVARIABLESTORAGE
SAVEOBJS = ADIGATORVARIABLESTORAGE.SAVE;
OVEROBJS = ADIGATORVARIABLESTORAGE.OVERMAP;
OVERLOCS = ADIGATOR.VARINFO.OVERMAP.IF;
SAVELOCS = ADIGATOR.VARINFO.SAVE.IF;

Remaps    = ADIGATOR.IFDATA(IfCount).BROS(BroCount).REMAPS;
IfOvermap = ADIGATOR.IFDATA(IfCount).OVERMAP;

if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.FORINFO.INNERLOC
  PriorDepStruc = ADIGATOR.IFDATA(IfCount).PRIORDEP;
end

for RemapLoc = Remaps
  if RemapLoc(2)
    % Remapping an object created within this conditional branch
    UnderObj = SAVEOBJS{RemapLoc(1)};
    OverObj  = OVEROBJS{OVERLOCS(RemapLoc(2))};
  else
    % Remapping an object created prior to conditional block
    PriorLoc = IfOvermap(RemapLoc(1)).Prior;
    if ADIGATOR.OPTIONS.UNROLL && ADIGATOR.FORINFO.INNERLOC && ...
        any(PriorDepStruc.Locs == PriorLoc)
      UnderObj = PriorDepStruc.Vars{PriorDepStruc.Locs == PriorLoc};
    else
      UnderObj = SAVEOBJS{SAVELOCS(PriorLoc,1)};
    end
    OverObj  = OVEROBJS{IfOvermap(RemapLoc(1)).Outer(1)};
  end
  cadaPrintReMap(UnderObj,OverObj,UnderObj.id);
end
end