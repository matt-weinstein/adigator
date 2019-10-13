function [ForCount, IfCount] = adigatorPrintTempFiles(Ffid,Tfid,FlowInfo,...
  DerNumber,ForCount,FunStrChecks)
% [ForCount, IfCount] = adigatorPrintTempFiles(Ffid,Tfid,FlowInfo,...
%   DerNumber,ForCount,FunStrChecks)
% This function is used to make blocks of temporary functions, which
% contain lines of user code interspersed with statements which call
% adigatorVarAnalyzer in order to read and print out the derivative function 
% properly.
% This routine is called from adigator.m and calls no routines which are not
% sub-routines of adigatorPrintTempFiles.m itself.
%
% Inputs:
% Ffid - user function file ID (reading from)
% Tfid - intermediate function file ID (writing to)
% FlowInfo - data structure containing file locations of all flow control
% DerNumber - integer identifying what order derivative file this is
% ForCount  - if non-zero, function is a not the primary function
% FunStrChecks - cell array containing strings to regexp to determine if a
%                function call exists on any given line
%
% Outputs:
% ForCount - number of for loops in the function
% IfCount  - number of if statements in the function
% 
%
% Copyright 2011-214 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

% ----------------------------------------------------------------------- %
%                                MAIN                                     %
% ----------------------------------------------------------------------- %

IfCount  = 0;
indent   = [];
MainFlag = ~ForCount;
ForIfFlag  = [0 0 0];
StartLocation = FlowInfo.StartLocation;
if ~isempty(FlowInfo.Children)
  % Print the first block of code.
  EndLocation = FlowInfo.Children(1).StartLocation;
  PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
  % Print the children blocks of code.
  [StartLocation,ForCount,IfCount,~] = FlowPrint(Ffid,FlowInfo.Children,...
    Tfid,indent,ForCount,IfCount,DerNumber,MainFlag,FunStrChecks,ForIfFlag);
end
EndLocation = FlowInfo.EndLocation;
% Print the last block of code.
PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
return
end

function [EndLocation,ForCount,IfCount,ForIfFlag] = FlowPrint(Ffid,FlowInfo,Tfid,...
  indent,ForCount,IfCount,DerNumber,MainFlag,FunStrChecks,ForIfFlag)
global ADIGATOR
% Function Recursively calls itself to get through the FlowInfo structure
% and print out all of the temporary files.
FlowSize = length(FlowInfo);
LastBro  = 0;
for FlowCount = 1:FlowSize
  if strcmp(FlowInfo(FlowCount).Type,'if')
    % ------------------------------------------------------------------- %
    %                         IF STATEMENTS - START                       %
    % ------------------------------------------------------------------- %
    IfCount       = IfCount+1;
    ForIfFlag(2)  = ForIfFlag(2)+1;
    BroCount = 1;
    CurIfCount = IfCount;
    if FlowCount ~= 1
      % Need to get the section between previous Flow Control Statment and
      % this one.
      StartLocation = EndLocation; % Starting at the End of the last section
      EndLocation = FlowInfo(FlowCount).StartLocation;
      PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
    end
    fprintf(Tfid,[indent,'%% ADiGator IF Statement #%1.0d: START\n'],IfCount);
    
    % --------------- Print Out the Conditional Variables --------------- %
    CondVarStr = getIfForStatement(Ffid,FlowInfo(FlowCount).StartLocation);
    CadaVarStr = 'cadaconditional1';
    if DerNumber == 1
      fprintf(Tfid,[indent,CadaVarStr,' = ',CondVarStr,';\n']);
      CondVarStr = FindDoinkers(CondVarStr);
      fprintf(Tfid,[indent,CadaVarStr,' = adigatorVarAnalyzer(''',CadaVarStr,...
        ' = ',CondVarStr,';'',',CadaVarStr,',''',CadaVarStr,''',0);\n']);
    end
    CondVarStrs = cell(1,FlowSize-FlowCount);
    CondVarStrs{1} = [CadaVarStr,','];
    % Get Any ElseIf variables
    for FlowCount2 = FlowCount+1:FlowSize
      if strcmp(FlowInfo(FlowCount2).Type,'elseif')
        BroCount = BroCount+1;
        if DerNumber == 1
          CondVarStr = getIfForStatement(Ffid,FlowInfo(FlowCount2).StartLocation);
          CadaVarStr = sprintf('cadaconditional%1.0d',BroCount);
          fprintf(Tfid,[indent,CadaVarStr,' = ',CondVarStr,';\n']);
          CondVarStr = FindDoinkers(CondVarStr);
          fprintf(Tfid,[indent,CadaVarStr,' = adigatorVarAnalyzer(''',CadaVarStr,...
            ' = ',CondVarStr,';'',',CadaVarStr,',''',CadaVarStr,''',0);\n']);
        end
        CondVarStrs{BroCount} = [CadaVarStr,','];
      elseif strcmp(FlowInfo(FlowCount2).Type,'else')
        BroCount = BroCount+1;
        CondVarStrs{BroCount} = '[],';
      else
        break
      end
    end
    LastBro = BroCount;
    % ---------------- Print Out the IF Initialize Statement ------------ %
    CondVarStrs = CondVarStrs(1:BroCount);
    CondVarStrs = cell2mat(CondVarStrs);
    if ADIGATOR.OPTIONS.UNROLL && (ForIfFlag(2) == 1 || ForIfFlag(3))
        fprintf(Tfid,[indent,'for adigatorIfPrint%1.0f = adigatorIfLooper(%1.0f)\n'],...
          IfCount,IfCount);
        fprintf(Tfid,[indent,'adigatorIfLooperi(adigatorIfPrint%1.0f,%1.0f);\n'],IfCount,IfCount);
    end
    fprintf(Tfid,[indent,'adigatorIfInitialize(%1.0f,',...
      CondVarStrs(1:end-1),');\n'],IfCount);
    if ADIGATOR.OPTIONS.UNROLL
      fprintf(Tfid,[indent,'[adigatorIfEvalStr, adigatorIfEvalVar] = ',...
        'adigatorIfIterStart(%1.0f,1);%%#ok<NASGU>\n'],IfCount);
      fprintf(Tfid,[indent,'if ~isempty(adigatorIfEvalStr)\nadigatorSetCellEvalFlag(1);',...
        ' cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);\nend\n']);
    else
      fprintf(Tfid,[indent,'adigatorIfIterStart(%1.0f,1);\n'],IfCount);
    end
    BroCount = 1;
  elseif strcmp(FlowInfo(FlowCount).Type,'elseif') || strcmp(FlowInfo(FlowCount).Type,'else')
    % ------------------------------------------------------------------- %
    %                  ELSEIF/ELSE STATEMENTS  - START                    %
    % ------------------------------------------------------------------- %
    BroCount = BroCount+1;
    % ------------------ Print Out the ELSEIF Statement ----------------- %
        fprintf(Tfid,[indent,'[adigatorIfEvalStr, adigatorIfEvalVar] = ',...
      'adigatorIfIterStart(%1.0f,%1.0f);%%#ok<NASGU>\n'],CurIfCount,BroCount);
    fprintf(Tfid,[indent,'if ~isempty(adigatorIfEvalStr)\nadigatorSetCellEvalFlag(1);',...
      ' cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);\nend\n']);
  elseif strcmp(FlowInfo(FlowCount).Type,'for')
    % ------------------------------------------------------------------- %
    %                      FOR STATEMENTS  - START                        %
    % ------------------------------------------------------------------- %
    ForCount      = ForCount+1;
    ForIfFlag(1)  = ForIfFlag(1)+1;
    if ForIfFlag(2); ForIfFlag(3) = ForIfFlag(3)+1; end
    CurForCount   = ForCount;
    if FlowCount ~= 1
      % Need to get the section between previous Flow Control Statment and
      % this one.
      StartLocation = EndLocation; % Starting at the End of the last section
      EndLocation = FlowInfo(FlowCount).StartLocation;
      PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
    end
    fprintf(Tfid,[indent,'%% ADiGator FOR Statement #%1.0d: START\n'],ForCount);
    % -------------------- Get the Loop Variable ------------------------ %
    [LoopStr,whileflag] = getIfForStatement(Ffid,FlowInfo(FlowCount).StartLocation);
    % See if we need to print it out or not
    if ADIGATOR.OPTIONS.UNROLL && whileflag
      error('Cannot unroll ''while'' loops - please set option to 0');
    end
    LoopVar    = sprintf('cadaforvar%1.0d',ForCount);
    LoopCount  = sprintf('cadaforcount%1.0d',ForCount);
    if whileflag
      if DerNumber == 1
        LoopVarStr = [LoopVar,' = ',LoopStr,';'];
        fprintf(Tfid,[indent,LoopVarStr,'\n']);
        fprintf(Tfid,[indent,LoopVar,' = adigatorVarAnalyzer(''',FindDoinkers(LoopVarStr),''',',...
          LoopVar,',''',LoopVar,''',0);\n']);
      else
        LoopVar = LoopStr;
      end
    else
      EqLoc = strfind(LoopStr,'=');
      if ~isempty(EqLoc)
        LoopStrLHS = strtrim(LoopStr(1:EqLoc(1)-1));
        LoopStrRHS = strtrim(LoopStr(EqLoc(1)+1:end));
        if DerNumber == 1
          LoopVarStr = [LoopVar,' = ',LoopStrRHS,';'];
          fprintf(Tfid,[indent,LoopVarStr,'\n']);
          fprintf(Tfid,[indent,LoopVar,' = adigatorVarAnalyzer(''',FindDoinkers(LoopVarStr),''',',...
            LoopVar,',''',LoopVar,''',0);\n']);
        else
          LoopVar = LoopStrRHS;
          LoopCount = LoopStrLHS;
        end
      else
        errlink = GenErrorLink(Ffid,FlowInfo(FlowCount).StartLocation(1));
        error(['???Unable to parse ',LoopStr,': No Equal Sign at: ',errlink])
      end
    end
    
    % --------------- Print the FOR Initialize Statement ---------------- %
    if ForIfFlag(1) == 1 && ~ADIGATOR.OPTIONS.UNROLL
      % Main function on an outer loop
      fprintf(Tfid,[indent,'[adigatorForVariable%1.0d, adigatorForEvalStr, adigatorForEvalVar]',...
        ' = adigatorForInitialize(%1.0d,',LoopVar,',%1.0f);%%#ok<NASGU>\n'],ForCount,ForCount,whileflag);
      fprintf(Tfid,[indent,'if ~isempty(adigatorForEvalStr)\n']);
      fprintf(Tfid,[indent,'    adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);\n']);
      fprintf(Tfid,[indent,'end\n']);
    else
      fprintf(Tfid,[indent,'adigatorForVariable%1.0d = adigatorForInitialize(%1.0d,',...
        LoopVar,',%1.0f);\n'],ForCount,ForCount,whileflag);
    end
    fprintf(Tfid,[indent,'for adigatorForVariable%1.0di = adigatorForVariable%1.0d\n'],...
      ForCount,ForCount);
    fprintf(Tfid,[indent,LoopCount,' = adigatorForIterStart(%1.0d,',...
      'adigatorForVariable%1.0di);\n'],ForCount,ForCount);
    if DerNumber == 1 && ~whileflag
      fprintf(Tfid,[indent,LoopStrLHS,' = ',LoopVar,'(:,',LoopCount,');\n']);
      fprintf(Tfid,[indent,LoopStrLHS,' = adigatorVarAnalyzer(''',LoopStrLHS,...
        ' = ',LoopVar,'(:,',LoopCount,');'',',LoopStrLHS,',''',LoopStrLHS,''',0);\n']);
    end
  end
  
  % --------------------------------------------------------------------- %
  %                      CALCULATIONS WITHIN THE FLOW CONTROL             %
  % --------------------------------------------------------------------- %
  indent = [indent,'    ']; %#ok<AGROW>
  StartLocation = FlowInfo(FlowCount).StartLocation;
  if ~isempty(FlowInfo(FlowCount).Children)
    % Print out the block leading to first child.
    EndLocation = FlowInfo(FlowCount).Children(1).StartLocation;
    PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber);
    % Print the children blocks of code.
    [StartLocation,ForCount,IfCount,ForIfFlag] = ...
      FlowPrint(Ffid,FlowInfo(FlowCount).Children,Tfid,indent,...
      ForCount,IfCount,DerNumber,MainFlag,FunStrChecks,ForIfFlag);
  end
  EndLocation = FlowInfo(FlowCount).EndLocation;
  % Print the last block of code for this section.
  PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber)

  % --------------------------------------------------------------------- %
  %                      END OF THE FLOW CONTROL BLOCK                    %
  % --------------------------------------------------------------------- %
  indent(1:4) = [];  
  if strcmp(FlowInfo(FlowCount).Type,'for')
    if whileflag && DerNumber == 1
      fprintf(Tfid,[indent,LoopVarStr,'\n']);
      fprintf(Tfid,[indent,LoopVar,' = adigatorVarAnalyzer(''',FindDoinkers(LoopVarStr),''',',...
        LoopVar,',''',LoopVar,''',0);\n']);
    end
    fprintf(Tfid,[indent,'[adigatorForEvalStr, adigatorForEvalVar]',...
      '= adigatorForIterEnd(%1.0d,adigatorForVariable%1.0di);\n'],CurForCount,CurForCount);
    fprintf(Tfid,[indent,'if ~isempty(adigatorForEvalStr)\n']);
    fprintf(Tfid,[indent,'    adigatorSetCellEvalFlag(1); cellfun(@eval,adigatorForEvalStr); adigatorSetCellEvalFlag(0);\n']);
    fprintf(Tfid,[indent,'end\n']);
    fprintf(Tfid,[indent,'end\n']);
    fprintf(Tfid,[indent,'%% ADiGator FOR Statement #%1.0d: END\n'],CurForCount);
    ForIfFlag(1) = ForIfFlag(1)-1;
    if ForIfFlag(2); ForIfFlag(3) = ForIfFlag(3) - 1; end
  elseif BroCount == LastBro
    fprintf(Tfid,[indent,'[adigatorIfEvalStr, adigatorIfEvalVar] = ',...
      'adigatorIfIterEnd(%1.0f,%1.0f);%%#ok<NASGU>\n'],CurIfCount,BroCount);
    fprintf(Tfid,[indent,'if ~isempty(adigatorIfEvalStr)\nadigatorSetCellEvalFlag(1);',...
      ' cellfun(@eval,adigatorIfEvalStr); adigatorSetCellEvalFlag(0);\nend\n']);
    if ADIGATOR.OPTIONS.UNROLL && (ForIfFlag(2) == 1 || ForIfFlag(3))
      fprintf(Tfid,[indent,'end\n']);
    end
    fprintf(Tfid,[indent,'%% ADiGator IF Statement #%1.0d: END\n'],CurIfCount);
    ForIfFlag(2) = ForIfFlag(2)-1;
  else
    fprintf(Tfid,[indent,'adigatorIfIterEnd(%1.0f,%1.0f);\n'],CurIfCount,BroCount);
  end
end

return
end

function PrintTemp(Ffid,StartLocation,EndLocation,Tfid,indent,FunStrChecks,DerNumber)
% Read from file Ffid, print to temporary function.
% Temporary function is named by the EndLocation(3).
global ADIGATOR
VAstr = 'adigatorVarAnalyzer';
fseek(Ffid,StartLocation(4),-1);

% ----------------First Line If There are Multiple Evaluations------------%
FunStrFULL = fgets(Ffid);
MajorLineCount = StartLocation(1);
MinorLineStart = StartLocation(2)+1;

while MajorLineCount <= EndLocation(1) && ~isnumeric(FunStrFULL)
  FunStrFull = strtrim(FunStrFULL);
  if ~isempty(FunStrFull)
    if length(FunStrFull) > 3 && ~strcmp(FunStrFull(1),'%')
      commentloc = findcomments(FunStrFull);
      if commentloc > 0
        FunStrFullTemp = strtrim(FunStrFull(1:commentloc-1));
        multlines1  = strcmp(FunStrFullTemp(end-2:end),'...');
        if multlines1
          FunStrFull = FunStrFullTemp;
        end
      else
        multlines1  = strcmp(FunStrFull(end-2:end),'...');
      end
      osquarelocs = strfind(FunStrFull,'[');
      csquarelocs = strfind(FunStrFull,']');
      multlines2  = length(osquarelocs) > length(csquarelocs);
      while multlines1 || multlines2
        % Single command spanning multiple lines - look for a comment at
        % end of this line
        commentloc = strfind(FunStrFull,'%');
        if ~isempty(commentloc)
          % Just trash the comment - no good way of keeping it
          FunStrFull = strtrim(FunStrFull(1:commentloc(1)-1));
        end
        if multlines1 && length(FunStrFull) > 6 && strcmp(FunStrFull(1:7),'global ')
          % ... at end of line with global
           FunStrFull = [FunStrFull(1:end-3),' '];
        elseif multlines1
          % ... at end of line
          FunStrFull = FunStrFull(1:end-3);
        elseif strcmp(FunStrFull(end),',')
          % Building a Matrix and have ',' at end for some reason
          FunStrFull = [FunStrFull(1:end-1),';'];
        elseif ~strcmp(FunStrFull(end),';')
          % Building a Matrix and dont have the vertcat operator
          FunStrFull = [FunStrFull,';']; %#ok<AGROW>
        end
        FunStrFull = strtrim([FunStrFull,fgets(Ffid)]);
        MajorLineCount = MajorLineCount + 1;
        commentloc = findcomments(FunStrFull);
        if commentloc > 0
          FunStrFullTemp = strtrim(FunStrFull(1:commentloc-1));
          multlines1  = strcmp(FunStrFullTemp(end-2:end),'...');
          if multlines1
            FunStrFull = FunStrFullTemp;
          end
        else
          multlines1  = strcmp(FunStrFull(end-2:end),'...');
        end
        osquarelocs = strfind(FunStrFull,'[');
        csquarelocs = strfind(FunStrFull,']');
        multlines2  = length(osquarelocs) > length(csquarelocs);
      end
    end
    [FunStr,NUMFunStr] = adigatorSeperateFunLines(FunStrFull);
    if MajorLineCount == EndLocation(1)
      NUMFunStr = EndLocation(2)-1;
    end
    % -----------------Work on Function Lines--------------------------
    for adigatorFScount = MinorLineStart:NUMFunStr
      FunStri = FunStr{adigatorFScount};
      StrLength = length(FunStri);
      FunStri = strtrim(FindComments(FunStri));
      SlashLocs = strfind(FunStri,'\');
      if ~isempty(SlashLocs)
        FunStri = strtrim(FindSlashes(FunStri,SlashLocs));
      end
      EqualLoc = regexp(FunStri,'[^=><~]=[^=]')+1;
      if strcmp(FunStri(1),'%')
        % COMMENT
        FunStri = FindDoinkers(FunStri);
        fprintf(Tfid,[indent,VAstr,'(''',FunStri,''');\n']);
      elseif StrLength >  6 && (strcmp(FunStri(1:6),'error(') || strcmp(FunStri(1:6),'error '))
        % ERROR MESSAGE
        IfCount  = evalin('caller','IfCount');
        BroCount = evalin('caller','BroCount');
        FunStri  = FindDoinkers(FunStri);
        fprintf(Tfid,[indent,'adigatorError(%1.0d,%1.0d,''',FunStri,''');\n'],...
          IfCount,BroCount);
      elseif ~isempty(EqualLoc)
        % Some sort of Assignment/Calculation
        FunStri = CheckFunctionCall(FunStri,FunStrChecks,Tfid,indent,DerNumber,MajorLineCount,Ffid);
        if isempty(FunStri)
          continue
        end
               
        % ---Work on LHS of equal Sign---
        VarStr = strtrim(FunStri(1:EqualLoc-1));
        if strcmp(VarStr(1),'[')
          % Multiple Assignments
          VarStrings = SeperateOutputStrings(strtrim(VarStr(2:end-1)),0);
          NumOutVars = length(VarStrings);
        else
          % Single Assignment
          VarStrings = cell(1); VarStrings{1} = VarStr;
          NumOutVars = 1;
        end
        SubsFlags = zeros(NumOutVars,1);
        
        %svacount = ADIGATOR.SVACOUNT;
        for Vcount = 1:NumOutVars
          VarStr = VarStrings{Vcount};
          SubsLoc = regexp(VarStr,'[\.\(\{]','once');
          if ~isempty(SubsLoc)
            SubsFlags(Vcount)  = 1;
            xStr = VarStr(1:SubsLoc-1);
            VarStrings{Vcount} = xStr;
            fprintf(Tfid,[indent,'if ~exist(''',xStr,''',''var''); ',...
              xStr,' = cadastruct([],''',xStr,''',[],0); end\n']);
            % Check for subs-assignment of structure or cell initialization
            RHSstr = strtrim(FunStri(EqualLoc+1:end));
            if any(~cellfun(@isempty,regexp(RHSstr,{'\<struct\(', '\<cell\(', '\<\{'})))
              % We want to make a temporary variable here.
              newVar = 'adigatorTempStruct';
              tempStr = [newVar,' = ',RHSstr];
              fprintf(Tfid,[indent,tempStr,'\n']);
              fprintf(Tfid,[indent,newVar,' = ',VAstr,'(''',tempStr,''',',...
                newVar,',''',newVar,''',0);\n']);
              FunStri = [FunStri(1:EqualLoc),newVar,';'];
            end
          end
        end
        
        %if DerNumber == 1 || svacount == ADIGATOR.SVACOUNT
        % Print Out the Actual Calculation
        fprintf(Tfid,[indent,FunStri,'\n']);
        
        % Get the Strings for the variables that we are going to feed into
        % the Variable Analyzer
        if NumOutVars == 1
          LHSout = VarStrings{1};
          RHSout = sprintf([VarStrings{1},',''',VarStrings{1},...
            ''',%1.0f'],SubsFlags(1));
        else
          LHSout = cell(1,NumOutVars); RHSout = cell(1,NumOutVars);
          LHSout{1} = ['[',VarStrings{1},','];
          RHSout{1} = sprintf([VarStrings{1},',''',VarStrings{1},...
            ''',%1.0f,'],SubsFlags(1));
          for Vcount = 2:NumOutVars-1
            LHSout{Vcount} = [VarStrings{Vcount},','];
            RHSout{Vcount} = sprintf([VarStrings{Vcount},',''',VarStrings{Vcount},...
              ''',%1.0f,'],SubsFlags(Vcount));
          end
          LHSout{NumOutVars} = [VarStrings{NumOutVars},']'];
          RHSout{NumOutVars} = sprintf([VarStrings{NumOutVars},',''',VarStrings{NumOutVars},...
            ''',%1.0f'],SubsFlags(NumOutVars));
          LHSout = cell2mat(LHSout); RHSout = cell2mat(RHSout);
        end
        % Print Out the call to the Variable Analyzer
        FunStri = FindDoinkers(FunStri);
        fprintf(Tfid,[indent,LHSout,' = ',VAstr,'(''',FunStri,''',',RHSout,');\n']);
      elseif StrLength > 7 && strcmp(FunStri(1:8),'keyboard')
        % KEYBOARD
        ADIGATOR.OPTIONS.KEYBOARD = 1;
        fprintf(Tfid,'keyboard\n');
        fprintf(Tfid,[indent,VAstr,'(''keyboard'');\n']);
      elseif StrLength > 4 && strcmp(FunStri(1:5),'break')
        % BREAK
        IfCount  = evalin('caller','IfCount');
        BroCount = evalin('caller','BroCount');
        fprintf(Tfid,[indent,'adigatorBreakCont(''break'',%1.0d,%1.0d);\n'],...
          IfCount,BroCount);
      elseif StrLength > 7 && strcmp(FunStri(1:8),'continue')
        % CONTINUE
        IfCount  = evalin('caller','IfCount');
        BroCount = evalin('caller','BroCount');
        fprintf(Tfid,[indent,'adigatorBreakCont(''continue'',%1.0d,%1.0d);\n'],...
          IfCount,BroCount);
      
      elseif StrLength > 5 && strcmp(FunStri(1:6),'return')
        % RETURN
        fprintf(Tfid,[indent,VAstr,'(''return'');\n']);
      elseif StrLength > 6 && strcmp(FunStri(1:7),'global ')
        % GLOBAL ASSIGNMENT
        VarString = strtrim(FunStri(8:end)); 
        if strcmp(VarString(end),';'); VarString = VarString(1:end-1); end
        VarLocs = strfind(VarString,' ');
        if ~isempty(VarLocs)
          % Check for multiple spaces.
          for Lcount = 2:length(VarLocs)
            if VarLocs(Lcount) == VarLocs(Lcount-1)+1
              VarLocs(Lcount-1) = 0;
            end
          end
          VarLocs = nonzeros(VarLocs);
          InVarStrs = cell(1,length(VarLocs)+1);
          OutVarStrs = cell(1,length(VarLocs)+1);
          InVarStrs{1} = ['''',strtrim(VarString(1:VarLocs(1)-1)),''','];
          OutVarStrs{1} = [strtrim(VarString(1:VarLocs(1)-1)),','];
          for Lcount = 2:length(VarLocs)
            InVarStrs{Lcount} = ['''',strtrim(VarString(...
              VarLocs(Lcount-1)+1:VarLocs(Lcount)-1)),''','];
            OutVarStrs{Lcount} = [strtrim(VarString(...
              VarLocs(Lcount-1)+1:VarLocs(Lcount)-1)),','];
          end
          InVarStrs{end} = ['''',strtrim(VarString(VarLocs(end)+1:end)),''''];
          OutVarStrs{end} = strtrim(VarString(VarLocs(end)+1:end));
          InVarStrs = cell2mat(InVarStrs);
          OutVarStrs = ['[',cell2mat(OutVarStrs),']'];
        else
          InVarStrs = ['''',VarString,''''];
          OutVarStrs = VarString;
        end
        fprintf(Tfid,[indent,OutVarStrs,' = ',VAstr,'(''global'',',InVarStrs,');\n']);
      elseif StrLength > 4 && strcmp(FunStri(1:5),'load(')
        % CHECK FOR LOAD - DONT ALLOW IT.
        errlink = GenErrorLink(Ffid,MajorLineCount);
        error(['In order to load in variables, must assign them to a'...
          ' variable name. At ',FunStri,' ',errlink]);
      elseif StrLength > 5 && strcmp(FunStri(1:6),'pause(')
        % PAUSE
      else
        % Dont know what this statement is.
        errlink = GenErrorLink(Ffid,MajorLineCount);
        error(['Cannot process statement: ',FunStrFull,' at ',errlink])
      end
    end
  end
  MinorLineStart = 1;
  FunStrFULL = fgets(Ffid);
  MajorLineCount = MajorLineCount+1;
end

return
end

function FunStri = CheckFunctionCall(FunStri,FunStrChecks,Tfid,indent,DerNumber,MajorLineCount,Ffid)

[Start,End] = regexp(FunStri,FunStrChecks,'start','end');

FunLoc = ~cellfun(@isempty,Start);
if FunLoc(1)
  errlink = GenErrorLink(Ffid,MajorLineCount);
  error(['User program cannot call main user function from within its methods at:' errlink])
elseif sum(FunLoc) > 1
  errlink = GenErrorLink(Ffid,MajorLineCount);
  error(['ADiGator only permits a single (user defined) function call per line. ',...
    'Please rewrite line: ',FunStri,' at: ',errlink])
elseif any(FunLoc)
  FunLocNum = 1:length(FunStrChecks);
  FunLoc = FunLocNum(FunLoc);
  Start = strfind(FunStri,FunStrChecks{FunLoc}(3:end));
  End   = End{FunLoc};
  if length(Start) > 1
    errlink = GenErrorLink(Ffid,MajorLineCount);
    error(['ADiGator only permits a single (user defined) function call per line. ',...
    'Please rewrite line: ',FunStri,' at: ',errlink])
  end
  fprintf(Tfid,[indent,'%% Call to User Function ',...
    FunStrChecks{FunLoc}(3:end-1),' --- (FunID %1.0d)\n'],FunLoc);
  % -------------------- Work on Inputs to Function --------------------- %
  InputStart = End+1; 
  % Get closing bracket for function call to get the input string
  InputEnd   = adigatorFindMatchingParen(FunStri,End)-1;
  InputStr = FunStri(InputStart:InputEnd);
  InputStrs = SeperateOutputStrings(InputStr,1);
  NumInputs = length(InputStrs);
  InVarStrs = cell(1,NumInputs);
  
  printNewIO = 0;
  for Icount = 1:NumInputs
    VarStr  = sprintf('cadainput%1.0d_%1.0d',FunLoc,Icount);
    InVarStrs{Icount} = [VarStr,';'];
    if ~strcmp(VarStr,InputStrs{Icount}) 
      printNewIO = 1;
      AsgnStr = [VarStr,' = ',InputStrs{Icount},';'];
      fprintf(Tfid,[indent,AsgnStr,'\n']);
      AsgnStr = FindDoinkers(AsgnStr);
      fprintf(Tfid,[indent,VarStr,' = adigatorVarAnalyzer(''',AsgnStr,''',',...
        VarStr,',''',VarStr,''',0);\n']);
    end
  end
  InVarStrs = cell2mat(InVarStrs);
  fprintf(Tfid,[indent,'adigatorInputs = {',InVarStrs(1:end-1),'};\n']);
  End = InputEnd+1;
  
  % -------------------- Print Call to Function ------------------------- %
  fprintf(Tfid,[indent,'[adigatorFunInfo, adigatorOutputs] = ',...
    'adigatortempfunc%1.0d(adigatorFunInfo,adigatorInputs);\n'],FunLoc);
  
  % ----------------------- Work On Outputs ----------------------------- %
  if ~isempty(regexp(FunStri(1:Start-1),'=\s*\>','once')) && strcmp(FunStri(1),'[') &&...
      (length(FunStri) == End+1 || ~isempty(regexp(FunStri(End+1:end),'<\\s*;','once')))
    % --------------------- Multiple Outputs ---------------------------- %
    EqLoc = strfind(FunStri,'=');
    OutputStart = 2;
    OutputEnd   = strfind(FunStri(1:EqLoc-1),']')-1;
    OutputStrs  = SeperateOutputStrings(FunStri(OutputStart:OutputEnd),0);
    NumOutput   = length(OutputStrs);
    for Ocount = 1:NumOutput
      OutVarStr = sprintf('cadaoutput%1.0d_%1.0d',FunLoc,Ocount);
      fprintf(Tfid,[indent,OutVarStr,' = adigatorOutputs{%1.0d};\n'],Ocount);
      if printNewIO && ~strcmp(OutputStrs{Ocount},'~')
        AsgnStr   = [OutputStrs{Ocount},' = ',OutVarStr,';'];
        fprintf(Tfid,[indent,AsgnStr,'\n']);
        SubLoc = regexp(OutputStrs{Ocount},'[\.\(\{]','once');
        if ~isempty(SubLoc)
          SubsStr = '1';
          VarStr  = OutputStrs{Ocount}(1:SubLoc(1)-1);
          fprintf(Tfid,[indent,'if ~exist(''',VarStr,''',''var''); ',...
            VarStr,' = cadastruct([],''',VarStr,''',[],0); end\n']);
        else
          SubsStr = '0';
          VarStr  = OutputStrs{Ocount};
        end
        fprintf(Tfid,[indent,VarStr,' = adigatorVarAnalyzer(''',AsgnStr,''',',...
          VarStr,',''',VarStr,''',',SubsStr,');\n']);
      end
    end
    FunStri = [];
  else
    % ----------------------- Single Output ----------------------------- %
    OutVarStr = sprintf('cadaoutput%1.0d_1',FunLoc);
    fprintf(Tfid,[indent,OutVarStr,' = adigatorOutputs{1};\n']);
    if printNewIO
      FunStri = [FunStri(1:Start-1),OutVarStr,FunStri(End+1:end)];
    else
      FunStri = [];
    end
  end
  
end
end

function FunStri = FindDoinkers(FunStri)

DoinkLocs = strfind(FunStri,'''');
if ~isempty(DoinkLocs)
  for Dcount = 1:length(DoinkLocs)
    Dloc = DoinkLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'''',FunStri(Dloc:end)];
  end
end
return
end

function FunStri = FindSlashes(FunStri,SlashLocs)

if ~isempty(SlashLocs)
  for Dcount = 1:length(SlashLocs)
    Dloc = SlashLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'\',FunStri(Dloc:end)];
  end
end
return
end

function FunStri = FindComments(FunStri)

DoinkLocs = strfind(FunStri,'%');
if ~isempty(DoinkLocs)
  for Dcount = 1:length(DoinkLocs)
    Dloc = DoinkLocs(Dcount)+Dcount;
    FunStri = [FunStri(1:Dloc-1),'%',FunStri(Dloc:end)];
  end
end
return
end

function VarStrings = SeperateOutputStrings(VarStr,IOflag)
% Just seperates output variable strings
% IOflag = 1 implies input
SpaceLocs = zeros(1,length(VarStr));
if ~IOflag
  % If outputs, treat spaces as seperators
  SpaceLocs(isspace(VarStr)) = 1;
end
CommaLocs = strfind(VarStr,',');
SpaceLocs(CommaLocs) = -1;
CharLocs  = 1:length(VarStr);

% Have a vector of locations where a zero or comma is.
ParenLoc1 = strfind(VarStr,'(');
if ~isempty(ParenLoc1)
  % Remove any entries of my zero/comma vector that are in
  % between parenthesis
  for Pcount = 1:length(ParenLoc1)
    CloseLoc = adigatorFindMatchingParen(VarStr,ParenLoc1(Pcount));
    SpaceLocs(CharLocs > ParenLoc1(Pcount) & CharLocs<CloseLoc) = 0;
  end
end
CurlyLoc1 = strfind(VarStr,'{');
if ~isempty(CurlyLoc1)
  % Remove any entries of my zero/comma vector that are in
  % between curlies
  for Ccount = 1:length(CurlyLoc1)
    CloseLoc = adigatorFindMatchingParen(VarStr,CurlyLoc1(Ccount));
    SpaceLocs(CharLocs > CurlyLoc1(Ccount) & CharLocs<CloseLoc) = 0;
  end
end
if IOflag
  % This is an input
  SquareLoc1 = strfind(VarStr,'[');
  if ~isempty(SquareLoc1)
    % Remove any entries of my zero/comma vector that are in
    % between square brackets
    for Ccount = 1:length(SquareLoc1)
      CloseLoc = adigatorFindMatchingParen(VarStr,SquareLoc1(Ccount));
      SpaceLocs(CharLocs > SquareLoc1(Ccount) & CharLocs<CloseLoc) = 0;
    end
  end
end
for Scount = 1:length(SpaceLocs)-1
  if SpaceLocs(Scount) && SpaceLocs(Scount+1)
    if SpaceLocs(Scount) < 1
      SpaceLocs(Scount+1) = 0;
      for S2count = Scount+2:length(SpaceLocs)-1
        if SpaceLocs(S2count) > 1
          SpaceLocs(S2count) = 1;
        else
          break
        end
      end
    else
      SpaceLocs(Scount) = 0;
    end
  end
end
SpaceLocs(end) = 0;
SpaceLocs1 = 1:length(VarStr);
SpaceLocs = SpaceLocs1(logical(SpaceLocs));

NumOutVars = length(SpaceLocs)+1;
VarLocs    = [0,SpaceLocs,length(VarStr)+1];
VarStrings = cell(NumOutVars,1);

for Vcount = 1:NumOutVars
  VarStrings{Vcount} = strtrim(VarStr(VarLocs(Vcount)+1:VarLocs(Vcount+1)-1));
end
end

function [Statement,whileflag] = getIfForStatement(Ffid,Location)
fseek(Ffid,Location(4),-1);
whileflag = 0;
StrFull = strtrim(fgets(Ffid));
while strcmp(StrFull(end-2:end),'...')
  StrFull = [StrFull(1:end-3),strtrim(fgets(Ffid))];
end

[MultStrs, ~] = adigatorSeperateFunLines(StrFull);
Statement     = MultStrs{Location(2)};
sLength       = length(Statement);
%Statement = FindDoinkers(Statement);
if sLength > 2 && strcmp(Statement(1:2),'if')
  Statement = strtrim(Statement(3:end));
elseif sLength > 6 && strcmp(Statement(1:6),'elseif')
  Statement = strtrim(Statement(7:end));
elseif sLength > 3 && strcmp(Statement(1:3),'for')
  Statement = strtrim(Statement(4:end));
elseif sLength > 5 && strcmp(Statement(1:5),'while')
  Statement = strtrim(Statement(6:end));
  whileflag = 1;
else
  errlink = GenErrorLink(Ffid,Location(1));
  error(['??? Unsure of what this statement is:',Statement,' at: ',errlink])
end
Statement = regexprep(Statement,'&&','&');
Statement = regexprep(Statement,'\|\|','|');
if strcmp(Statement(end),';')
  Statement = strtrim(Statement(1:end-1));
end
if strcmp(Statement(end),',')
  Statement = strtrim(Statement(1:end-1));
end
end

function errlink = GenErrorLink(Ffid,MajorLineCount)
filename     = fopen(Ffid);
filenameloc  = strfind(filename,filesep);
filenameDisp = filename(filenameloc(end)+1:end);
filename     = regexprep(filename,'\\','\\\\');
filenameDisp = regexprep(filenameDisp,'\\','\\\\');
errlink = sprintf(['User''s Function: <a href="matlab: opentoline(',filename,',%1.0d)">',filenameDisp,' line %1.0d</a>'],MajorLineCount,MajorLineCount);
end

function comment = findcomments(str)
% Look for a comment at end of string
comment = 0;
commentlocs = strfind(str,'%');
if ~isempty(commentlocs)
  % Make sure no %'s are used to build a string
  leftstrlocs = regexp(str,'[=''\(,]\s*''');
  rightstrlocs = regexp(str,'''\s*[''\),]');
  % if these lengths arent equal its possible that they occur after the
  % comment
  if length(leftstrlocs) < length(rightstrlocs)
    rightstrlocs = rightstrlocs(1:length(leftstrlocs));
  elseif length(leftstrlocs) > length(rightstrlocs)
    leftstrlocs = leftsrlocs(1:length(rightstrlocs));
  end
  if ~isempty(leftstrlocs)
    for C = commentlocs
      if ~any(C > leftstrlocs & C < rightstrlocs)
        comment = C;
      end
    end
  else
    comment = commentlocs(1);
  end
end
end
