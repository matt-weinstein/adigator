function [myLoopVar, outEvalStr, outEvalVar] = adigatorForInitialize(ForCount,UserLoopVar,whileflag)
% This transformation routine is called prior to the evaluation of any FOR
% loop in the intermediate program.
%
% Inputs:
%   ForCount - integer identifying the FOR loop
%   UserLoopVar - the variable in the user's program to be looped over
%   whileflag - binary variable (true means this is actually a WHILE loop
%               UserLoopVar is the arguement of the WHILE loop)
%
% Outputs:
%   myLoopVar - the variable to be looped over in the intermediate program
%   outEvalStr - cell arary of strings to be evaluated on the output in 
%                order to modify the workspace
%   outEvalVar - cell array containing variables to be placed into the 
%                workspace 
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE
maxWhileIters = ADIGATOR.OPTIONS.MAXWHILEITER;
if ~ADIGATOR.RUNFLAG
  % --------------------------------------------------------------------- %
  %                            Empty Run                                  %
  % --------------------------------------------------------------------- %
  ADIGATORFORDATA(ForCount).COUNTNAME = sprintf('cadaforcount%1.0d',ForCount);
  
  ADIGATORFORDATA(ForCount).START     = ADIGATOR.VARINFO.COUNT;
  ADIGATORFORDATA(ForCount).PARENTLOC = ADIGATOR.FORINFO.INNERLOC;
  ADIGATORFORDATA(ForCount).WHILEFLAG = whileflag;

  if ADIGATOR.FORINFO.FLAG
    ADIGATORFORDATA(ADIGATOR.FORINFO.INNERLOC).CHILDLOCS(end+1,1) = ForCount;
  elseif ~ADIGATOR.FORINFO.OUTERLOC
    ADIGATOR.FORINFO.OUTERLOC = ForCount;
  end
  ADIGATOR.FORINFO.INNERLOC = ForCount;
  ADIGATORFORDATA(ForCount).CHILDLOCS = zeros(0,1);
  if ~ADIGATOR.OPTIONS.UNROLL
    ADIGATOR.FORINFO.FLAG = 1;
  end
  outEvalStr = [];
  outEvalVar = [];
  myLoopVar = 0;
  
elseif ADIGATOR.OPTIONS.UNROLL
  % --------------------------------------------------------------------- %
  %            Unrolling Loop - Overmapping or Printing Eval              %
  % --------------------------------------------------------------------- %
  if ADIGATOR.EMPTYFLAG
    % Dont run this loop
    myLoopVar = [];
  elseif isa(UserLoopVar,'cada')
    ForLength = UserLoopVar.func.size(2);
    if isinf(ForLength)
      error('Cannot loop over vectorized dimension')
    end
    myLoopVar = 1:ForLength;
  else
    myLoopVar = UserLoopVar;
    ForLength = length(myLoopVar);
  end
  % ---------------------- Set For Lengths ---------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;

  if ParentLoc
    ParentIter = ADIGATORFORDATA(ParentLoc).COUNT.ITERATION;
    if ParentIter == 1
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
    else
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(ParentIter) = ForLength;
    end
  else
    ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
  end
  
  if isempty(myLoopVar)
    % Dont need to run loop
    ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).END+1;
    ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  else
    ADIGATOR.FORINFO.INNERLOC = ForCount;
    if ~ParentLoc
      ADIGATOR.FORINFO.OUTERLOC = ForCount;
    end
  end
elseif ADIGATOR.RUNFLAG == 1
  % --------------------------------------------------------------------- %
  %                           OverMap Run                                 %
  % --------------------------------------------------------------------- %
  if ADIGATOR.EMPTYFLAG
    % Dont run this loop
    ForLength = 0;
    myLoopVar = [];
  elseif whileflag
    ForLength = maxWhileIters;
    myLoopVar = 1:maxWhileIters;
  elseif isa(UserLoopVar,'cada')
    ForLength = UserLoopVar.func.size(2);
    if isinf(ForLength)
      error('Cannot loop over vectorized dimension')
    end
    myLoopVar = 1:ForLength;
  else
    myLoopVar = UserLoopVar;
    ForLength = length(myLoopVar);
  end
  % ---------------------- Set For Lengths ---------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;
  if ParentLoc
    ParentIter = ADIGATORFORDATA(ParentLoc).COUNT.ITERATION;
    if ParentIter == 1
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
    else
      ADIGATORFORDATA(ForCount).FOR(1).LENGTHS(ParentIter) = ForLength;
    end
  else
    ADIGATORFORDATA(ForCount).FOR(1).LENGTHS = ForLength;
  end
  
  if isempty(myLoopVar)
    % Dont ned to run this loop. Set exit sequence
    ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).END+1;
    ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  else
    % -------------- Overmap Variables Coming into Loop ----------------- %
    ADIGATOR.FORINFO.INNERLOC = ForCount;
    if ~ParentLoc
      ADIGATOR.FORINFO.OUTERLOC = ForCount;
    end
    if whileflag
      % Get the inputs for WHILE loop
      whilelocs = zeros(1,length(ADIGATORFORDATA(ForCount).PREVOVERMAP));
      whilevars = cell(1,length(ADIGATORFORDATA(ForCount).PREVOVERMAP));
    end
    count = 0;
    for VarCount = ADIGATORFORDATA(ForCount).PREVOVERMAP
      count= count+1;
      if VarCount < ADIGATOR.VARINFO.COUNT
        % This variable gets created prior to the loop
        SaveLoc = ADIGATOR.VARINFO.SAVE.FOR(VarCount,1);
        PrevVar = ADIGATORVARIABLESTORAGE.SAVE{SaveLoc};
        OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(VarCount,2);
        OverVar = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
      else
        % This variable is getting subsasgn'd to without pre-allocating.
        PrevVar = [];
        OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(VarCount,1);
        OverVar = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
      end
      if ~isa(OverVar,'cada') && ~isa(OverVar,'cadastruct') && isempty(OverVar)
        OverVar = PrevVar;
      else
        OverVar = cadaUnionVars(PrevVar,OverVar);
      end
      ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc} = OverVar;
      if whileflag
        % Collect while inputs
        whilelocs(count) = OverLoc;
        whilevars{count} = OverVar;
      end
    end
    if whileflag
      ADIGATORFORDATA(ForCount).WHILEINPUTS.LOCS = whilelocs;
      ADIGATORFORDATA(ForCount).WHILEINPUTS.VARS = whilevars;
    end
    % Checks will be done in adigatorForIterEnd to determine when all
    % inputs structures have become static..
    ADIGATOR.FORINFO.FLAG = 1;
  end
  outEvalStr = [];
  outEvalVar = [];
  
elseif ADIGATORFORDATA(ForCount).MAXLENGTH == 0
  % --------------------------------------------------------------------- %
  %                  Printing Run - Loop Never Runs                       %
  % --------------------------------------------------------------------- %
  % Dont ned to run this loop. Set exit sequence
  ADIGATOR.VARINFO.COUNT = ADIGATORFORDATA(ForCount).END+1;
  ADIGATOR.PREOPCOUNT    = ADIGATOR.VARINFO.COUNT;
  myLoopVar = [];
  outEvalVar = [];
  outEvalStr = [];
else
  % --------------------------------------------------------------------- %
  %                           Printing Run                                %
  % --------------------------------------------------------------------- %
  ParentLoc = ADIGATORFORDATA(ForCount).PARENTLOC;
  ADIGATOR.FORINFO.INNERLOC = ForCount;
  ADIGATOR.FORINFO.FLAG     = 1;
  fid    = ADIGATOR.PRINT.FID;
  indent = ADIGATOR.PRINT.INDENT;
  if ~ParentLoc
    % Is an Outer Loop
    ADIGATOR.FORINFO.OUTERLOC = ForCount;
    % -------------- Remap Variables Coming Into Loop ------------------- %
    nOutEval = length(ADIGATORFORDATA(ForCount).PREVOVERMAP);
    outEvalStr = cell(nOutEval,1);
    outEvalVar = cell(nOutEval,1);
    for Vcount = 1:nOutEval
      % Do Re-Map
      VarCount = ADIGATORFORDATA(ForCount).PREVOVERMAP(Vcount);
      if VarCount < ADIGATOR.VARINFO.COUNT
        % This variable is assigned prior to loop
        OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(VarCount,2);
        PrevVar  = ADIGATORVARIABLESTORAGE.SAVE{ADIGATOR.VARINFO.SAVE.FOR(VarCount,1)};
      else
        % This is an un-allocated subsasgn
        OverLoc = ADIGATOR.VARINFO.OVERMAP.FOR(VarCount,1);
        PrevVar = [];
      end
      NameLoc  = ADIGATOR.VARINFO.NAMELOCS(VarCount,1);
      VarStr   = ADIGATOR.VARINFO.NAMES{NameLoc};
      
      
      
      if OverLoc
        OverVar  = ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc};
        OverVar  = cadaPrintReMap(PrevVar,OverVar,VarCount);
        ADIGATORVARIABLESTORAGE.OVERMAP{OverLoc} = OverVar;
        outEvalVar{Vcount} = OverVar;
      else
        outEvalVar{Vcount} = PrevVar;
      end
      outEvalStr{Vcount} = sprintf([VarStr,' = adigatorForEvalVar{%1.0d};'],Vcount);
      dotLocs  = strfind(VarStr,'.');
      if ~isempty(dotLocs)
        StrucStr = VarStr(1:dotLocs(1)-1);
        strucEvalStr = ['if ~exist(''',StrucStr,''',''var''); ',StrucStr,...
          ' = cadastruct([],''',StrucStr,''',[],0); end;'];
        outEvalStr{Vcount} = [strucEvalStr,outEvalStr{Vcount}];
      end
      % Check for case when conditional set needs to return a saved
      % variable.
      if ~isempty(ADIGATOR.VARINFO.SAVE.IF)
        IfSaveLoc = ADIGATOR.VARINFO.SAVE.IF(VarCount,1);
        if IfSaveLoc
          ADIGATORVARIABLESTORAGE.SAVE{IfSaveLoc} = outEvalVar{Vcount};
        end
      end
    end
    
    % ---------------------- Print the Loop ----------------------------- %
    if ~ADIGATOR.EMPTYFLAG
      if whileflag
        fprintf(fid,[indent,'while ',UserLoopVar.func.name,';\n']);
      else
        fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
          ' = 1:%1.0d\n'],ADIGATORFORDATA(ForCount).MAXLENGTH);
      end
    end
  else
    % Is an Inner Loop
    % -------Print Out any Organizational Function Data Reshape/Refs----- %
    OrgFuncs = {'SUBSREF','SUBSASGN','SPARSE','NONZEROS','HORZCAT',...
      'VERTCAT','TRANSPOSE','REPMAT','RESHAPE','SIZE','STRUCTREF','STRUCTASGN'};
    for OFcount = 1:length(OrgFuncs)
      DataStructure = ADIGATORFORDATA(ForCount).(OrgFuncs{OFcount});
      % do .INDICES fields first
      if strcmp(OrgFuncs{OFcount},'HORZCAT') ||...
          strcmp(OrgFuncs{OFcount},'VERTCAT') || ...
          strcmp(OrgFuncs{OFcount},'STRUCTREF') || ...
          strcmp(OrgFuncs{OFcount},'STRUCTASGN')
        % horzcat and vertcat have added dimension
        for OF2count = 1:length(DataStructure)
          StringArray = DataStructure(OF2count).INDICES;
          for Icount = 1:size(StringArray,1)
            for Jcount = 1:size(StringArray,2)
              LHSstr = StringArray{Icount,Jcount,1};
              RHSstr = StringArray{Icount,Jcount,2};
              if ~isempty(LHSstr) && ~isempty(RHSstr)
                fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
              end
            end
          end
        end
      elseif isfield(DataStructure,'INDICES')
        % will get the other ones
        for OF2count = 1:length(DataStructure)
          StringArray  = DataStructure(OF2count).INDICES;
          SubsAsgnFlag = size(StringArray,2)>3;
          for Icount = 1:size(StringArray,1)
            LHSstr = StringArray{Icount,1};
            RHSstr = StringArray{Icount,2};
            if ~isempty(LHSstr) && ~isempty(RHSstr)
              fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
            end
            if SubsAsgnFlag % checks for subsasgn, he has extra
              LHSstr = StringArray{Icount,4};
              RHSstr = StringArray{Icount,5};
              if ~isempty(LHSstr) && ~isempty(RHSstr)
                fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
              end
            end
          end
        end
      end
      if isfield(DataStructure,'SIZES')
        % do .SIZES field next - all are straightforward
        for OF2count = 1:length(DataStructure)
          StringArray = DataStructure(OF2count).SIZES;
          if iscell(StringArray)
            for Icount = 1:size(StringArray,1)
              LHSstr = StringArray{Icount,1};
              RHSstr = StringArray{Icount,2};
              if ~isempty(LHSstr) && ~isempty(RHSstr)
                fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
              end
            end
          end
        end
      end
    end

    % print out any FOR loop changing size stuff
    for Fcount = 1:length(ADIGATORFORDATA(ForCount).FOR)
      LHSstr = ADIGATORFORDATA(ForCount).FOR(Fcount).LENGTHS{1};
      RHSstr = ADIGATORFORDATA(ForCount).FOR(Fcount).LENGTHS{2};
      if ~isempty(LHSstr) && ~isempty(RHSstr)
        fprintf(fid,[indent,LHSstr,' = ',RHSstr,';\n']);
      end
    end
    
    % -------------------- Print the Loop ------------------------------- %
    if ~ADIGATOR.EMPTYFLAG
      if ADIGATOR.DERNUMBER == 1
        LoopVarStr = ADIGATORFORDATA(ForCount).FOR(1).LENGTHS{1};
        if ~isempty(LoopVarStr)
          % Size of loop is dependent upon an outer loop
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:',LoopVarStr,'\n']);
        else
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:%1.0d\n'],ADIGATORFORDATA(ForCount).MAXLENGTH);
        end
      else
        if isa(UserLoopVar,'cada')
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:',UserLoopVar.func.name,'\n']);
        else
          fprintf(fid,[indent,'for ',ADIGATORFORDATA(ForCount).COUNTNAME,...
            ' = 1:%1.0d\n'],length(UserLoopVar));
        end
      end
    end
    outEvalStr = [];
    outEvalVar = [];
  end
  ADIGATOR.PRINT.INDENT = [indent,'    '];
  myLoopVar = 1;
end

if ADIGATOR.OPTIONS.UNROLL
  ADIGATOR.FORINFO.FLAG = 0;
end