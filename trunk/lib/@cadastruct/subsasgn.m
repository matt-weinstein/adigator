function x = subsasgn(x,s,b)
% CADASTRUCT overloaded SUBSASGN
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0

global ADIGATOR
nsubs = length(s);
y  = x.val;
x0 = y;
newfieldflag = 0;
if nsubs == 1 && strcmp(s.type,'()') && (isa(b,'cada') || isnumeric(b)) && isnumeric(x.val) && isempty(x.val)
  x = [];
  x = subsasgn(x,s,b);
  return
end
if isstruct(b) || iscell(b)
  error(['Cannot parse RHS of assignment - ',...
    'please instantiate the RHS cell/structure by directly assigning to a variable.']);
end
if ~x.arrayflag
  % --------------------------------------------------------------------- %
  %                           Scalar Structure                            %
  % --------------------------------------------------------------------- %
  % Parse through first
  refstr = cell(1,nsubs);
  for subcount = 1:nsubs
    i = s(subcount);
    if strcmp(i.type,'.')
      if ~newfieldflag && isfield(y,i.subs)
        y = subsref(y,i);
      elseif ~newfieldflag
        newfieldflag = subcount;
        y = [];
      else
        y = [];
      end
      refstr{subcount} = ['.',i.subs];
    elseif strcmp(i.type,'()')
      dotcount = subcount - 1;
      if isstruct(y) || nsubs > subcount || ...
          (newfieldflag && (isstruct(b) || iscell(b) || isa(b,'cadastruct')))
        % Need to convert to a structure array.
        yname = [x.name,cell2mat(refstr(1:dotcount))];
        if ~isempty(y)
          count = structconvert(y,yname,0);
          if ~count
            count = x.id;
          end
        else
          count = x.id;
        end
        yname = [x.name,cell2mat(refstr(1:dotcount))];
        y = cadastruct(y,yname,count,1);
      end
      break
    else
      dotcount = subcount - 1;
      break
    end
    if subcount == nsubs
      dotcount = nsubs;
      break
    end
    if isa(y,'cada') && subcount < nsubs && ~(subcount == nsubs-1 &&...
        strcmp(s(nsubs).type,'()') && (isnumeric(b) || isa(b,'cada'))) 
      % Empty cada object, remove it.
      y = [];
      x.val = subsasgn(x.val,s(1:subcount),[]);
      x0 = x.val;
    elseif isa(y,'cada') || isa(y,'cadastruct')
      dotcount = subcount;
      break
    end
  end

  if dotcount < nsubs
    % Cant do all of the work here, need to call structure array subsasgn
    % or cada subsasgn
    if isnumeric(y) && isempty(y)
      snew = s(dotcount+1:end);
      if length(snew) > 1 || ~strcmp(snew(1).type,'()') || isa(b,'cadastruct')
        yname = [x.name,cell2mat(refstr(1:dotcount))];
        y = cadastruct([],yname,[],1);
      end
      b = subsasgn(y,snew,b);
    else
      b = subsasgn(y,s(dotcount+1:end),b);
    end
  end

  if dotcount > 0
    
    % Call VarAnalyzer on b with the new name.
    bStr = [x.name,cell2mat(refstr(1:dotcount))];
    if ~ADIGATOR.CELLEVALFLAG
      if isa(b,'cadastruct') && dotcount == nsubs && ADIGATOR.RUNFLAG == 2
        FunString = [bStr,' = ',b.name,';'];
      elseif isnumeric(b)
        % b is numeric input
        [bMrow,bNcol] = size(b);
        btemp.func = struct('name',[],'size',[bMrow,bNcol],'zerolocs',[],...
          'value',b);
        if ADIGATOR.PRINT.FLAG
          if bMrow*bNcol < 10 && all(floor(b(:)) == b(:))
            btemp.func.name = mat2str(b);
          else
            btemp.func.name = cadamatprint(b);
          end
        end
        btemp.func.size = [bMrow, bNcol];
        NUMvod  =   ADIGATOR.NVAROFDIFF;
        btemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
        b = cada([],btemp.func,btemp.deriv);
        %FunString = [bStr,' = ',b.func.name,';'];
        FunString = '';
      elseif ischar(b)
        FunString = [bStr,' = ''',b,''';'];
      else
        FunString = '';
      end
      b = adigatorVarAnalyzer(FunString,b,bStr,dotcount < nsubs);
    end
    if isa(b,'cadastuct') && ~b.arrayflag
      b = b.val;
    end
    
    % Assign b to x.
    if isa(b,'cada') || isa(b,'cadastruct')
      % MATLAB will call subsasgn recursively instead of just doing the
      % structure/array assignment if use subsasgn(a,s,b) and b is
      % overloaded and a is cell/structure - this will NOT call overloaded
      % subsasgn if use a.field = b or a{i} = b which is super strange -
      % this is a kluge around it.
      if dotcount > 1
        if newfieldflag > 0 && newfieldflag < dotcount
          y = [];
        else
          y = subsref(x0,s(1:dotcount-1));
        end
      else
        y = x.val;
      end
      if isa(y,'cada')
        y = [];
      end
      y.(s(dotcount).subs) = b;
      if dotcount > 1
        x0 = subsasgn(x0,s(1:nsubs-1),y);
      else
        x0 = y;
      end
    else
      x0 = subsasgn(x0,s(1:dotcount),b);
    end
    x.val = x0;
  else
    % Had to convert to a structure array
    x = b;
  end
else
  % --------------------------------------------------------------------- %
  %                         Structure/Cell Array                          %
  % --------------------------------------------------------------------- %
  NUMvod  =   ADIGATOR.NVAROFDIFF;
  fid     =   ADIGATOR.PRINT.FID;
  PFLAG   =   ADIGATOR.PRINT.FLAG;
  indent  =   ADIGATOR.PRINT.INDENT;
  
  if ~isa(b,'cada') && ~isa(b,'cadastruct') && isnumeric(b) && any(size(b)) > 0
    % b is numeric input - give it its own ID so that we can track it
    % better..
    [bMrow,bNcol] = size(b);
    btemp.func = struct('name',[],'size',[bMrow,bNcol],'zerolocs',[],...
      'value',b);
    bid = ADIGATOR.VARINFO.COUNT;
    ADIGATOR.VARINFO.LASTOCC(bid,1) = ADIGATOR.VARINFO.COUNT;
    if PFLAG
      bname = cadafuncname();
      btemp.func.name = bname;
      if bMrow*bNcol < 10 && all(floor(b(:)) == b(:))
        fprintf(fid,[indent,bname,' = ',mat2str(b),';\n']);
      else
        fprintf(fid,[indent,bname,' = ',cadamatprint(b),';\n']);
      end
    end
    btemp.func.size = [bMrow, bNcol];
    btemp.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    b = cada(bid,btemp.func,btemp.deriv);
    ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  end
  
  if ADIGATOR.FORINFO.FLAG
    IncreaseForAsgnCount();
    if ADIGATOR.RUNFLAG == 2
      [y,returnflag] = ForSubsAsgn(x,s,b);
      if returnflag
        x = y;
        return
      end
      y = x.val;
    end
  end
  
  
  if (isa(b,'cada') || isnumeric(b)) && strcmp(s(end).type,'()')
    refasgnflag = 1;
    nsubs = nsubs-1;
  else
    refasgnflag = 0;
  end
  % Make an initial sweep through - build reference string, replace cada
  % indices with numeric, check for empty fields which are cada objects and
  % replace with empty arrays
  snumeric = s;
  cadaflag = 0;
  refstr = cell(1,nsubs);
  for subcount = 1:nsubs
    i = s(subcount);
    if isa(y,'cada') && prod(y.func.size) == 0
      y = [];
    end
    switch i.type
      case '.'
        refstr{subcount} = ['.',i.subs];
        if ~newfieldflag && isfield(y,i.subs)
          if numel(y) > 1 && ADIGATOR.EMPTYFLAG
            y = y(1);
            if subcount > 1
              x.val = subsasgn(x.val,s(1:subcount-1),y);
            else
              x.val = y;
            end
          end
          y = subsref(y,i);
        else
          if ~newfieldflag
             newfieldflag = subcount;
          end
          y = [];
        end
      case '()'
        [refstri,i.subs,cadaflag2] = parseIndex(i.subs);
        cadaflag = or(cadaflag,cadaflag2);
        snumeric(subcount) = i;
        refstr{subcount} = ['(',refstri,')'];
        if ~newfieldflag
          try
            y = subsref(y,i);
          catch
            y(:) = [];
            newfieldflag = subcount;
          end
        else
          y(:) = [];
        end
      case '{}'
        [refstri,i.subs,cadaflag2] = parseIndex(i.subs);
        cadaflag = or(cadaflag,cadaflag2);
        snumeric(subcount) = i;
        refstr{subcount} = ['{',refstri,'}'];
        if ~newfieldflag
          try
            y = subsref(y,i);
          catch
            y = [];
            newfieldflag = subcount;
          end
        end
      otherwise
        error('?? invalid assignment type')
    end
  end
  

  yprior = y;
  
  if refasgnflag
    if isnumeric(b) && isempty(b) && ~newfieldflag && (isstruct(y) || iscell(y))
      % Isn't actually a cada object, rather it is removing an element of a
      % structure/cell array
      nsubs = nsubs+1;
      i = s(nsubs);
      [refstri,i.subs,cadaflag2] = parseIndex(i.subs);
      cadaflag = or(cadaflag,cadaflag2);
      refstr{nsubs} = ['(',refstri,')'];
      snumeric(nsubs) = i;
      refasgnflag = 0;
      asgnstring = [x.name,cell2mat(refstr)];
    else
      % do assignment to b.
      asgnstring = [x.name,cell2mat(refstr)];

      if isa(y,'cada')
        y.id = x.id;
      end
      if PFLAG && ~ADIGATOR.EMPTYFLAG
        % We want CADA subsasgn to print out the assignment - we are going
        % to modify the VARINFO.NAMES corresponding to this variable so
        % that it gets the names properly.
        
        count = ADIGATOR.VARINFO.COUNT;
        nameloc = ADIGATOR.VARINFO.NAMELOCS(count,1);
        oldname = ADIGATOR.VARINFO.NAMES{nameloc};
        ADIGATOR.VARINFO.NAMES{nameloc} = asgnstring;
        % Rename y
        if isa(y,'cada')
          funcstr = cadafuncname();
          yderiv = y.deriv;
          for Vcount = 1:NUMvod
            if ~isempty(yderiv(Vcount).nzlocs)
              yderiv(Vcount).name = cadadername(funcstr,Vcount);
            end
          end
          y.deriv = yderiv;
          y.func.name = funcstr;
        end
        % Call subsasgn
        b = subsasgn(y,s(end),b);
        % Replace VARINFO.NAMES
        ADIGATOR.VARINFO.NAMES{nameloc} = oldname;
      else
        b = subsasgn(y,s(end),b);
      end
      % Remove the count for this.
      ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT-1;
    end
  else
    asgnstring = [x.name,cell2mat(refstr)];
  end
  
  % ------------------------ Do Assignment ------------------------------ %
  xoldid = x.id;
  xid  = ADIGATOR.VARINFO.COUNT;
  xnew.id   = xid;
  xnew.name = x.name;
  
  
  if isa(b,'cada')
    % MATLAB behaves very oddly when b is a cada and we use subsasgn(a,s,b)
    % and behaves differently if we use a.x = b - very strange stuff and a
    % major pain..
    
    if nsubs > 1
      if newfieldflag == nsubs-1 && strcmp(snumeric(nsubs-1).type,'()')
        if nsubs > 2
          y = subsref(x.val,snumeric(1:nsubs-2));
        else
          y = x.val;
        end
        if isstruct(y)
          fnames = fieldnames(y);
          structinit = cell(2,length(fnames));
          structinit(1,:) = fnames;
          y = struct(structinit{:});
        else
          y = [];
        end
      elseif newfieldflag < nsubs && newfieldflag > 0
        y = [];
      else
        y = subsref(x.val,snumeric(1:nsubs-1));
      end
    else
      y = x.val;
    end
    i = snumeric(nsubs);
    if isa(y,'cada')
      y = [];
    end
    switch i.type
      case '.'
        y(1).(i.subs) = b;
      case '{}'
        y{i.subs{:}} = b;
      case '()'
        y(i.subs{:}) = b;
    end
    if nsubs > 1
      if strcmp(snumeric(1).type,'()') && isnumeric(x.val) && isempty(x.val)
        x.val = struct(snumeric(2).subs,[]);
      elseif strcmp(snumeric(nsubs-1).type,'()') && isstruct(y) && newfieldflag
        % x(1:nsubs-1) is a structure - need to protect against dissimilar structure error
        x.val = subsasgn(x.val,snumeric(1:nsubs),[]);
      end
      xnew.val = subsasgn(x.val,snumeric(1:nsubs-1),y);
    else
      xnew.val = y;
    end

    if ~refasgnflag
      % Only do this if there's no call to CADA subsasgn
      if PFLAG && ~ADIGATOR.EMPTYFLAG
        nameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
        oldname = ADIGATOR.VARINFO.NAMES{nameloc};
        ADIGATOR.VARINFO.NAMES{nameloc} = asgnstring;
        funasgnstr = cadafuncname(xid);
      end
      
      if ADIGATOR.RUNFLAG && ~ADIGATOR.EMPTYFLAG
        dbflag = cadaCheckForDerivs(b);
      else
        dbflag = 0;
      end
      
      if PFLAG && dbflag
        [~,DPFLAG] = cadafuncname(b.id);
      else
        DPFLAG = 0;
      end
      
      if DPFLAG && ~ADIGATOR.EMPTYFLAG
        bderiv = b.deriv;
        for Vcount = 1:NUMvod
          if ~isempty(bderiv(Vcount).nzlocs)
            derasgnstr = cadadername(funasgnstr,Vcount,xid);
            fprintf(fid,[indent,derasgnstr,' = ',bderiv(Vcount).name,';\n']);
          end
        end
      end
      
      if PFLAG && ~ADIGATOR.EMPTYFLAG
        fprintf(fid,[indent,funasgnstr,' = ',b.func.name,';\n']);
        ADIGATOR.VARINFO.NAMES{nameloc} = oldname;
      end
    end
    ADIGATOR.VARINFO.LASTOCC(b.id,1) = ADIGATOR.VARINFO.COUNT;
  elseif isa(b,'cadastruct')
    if ADIGATOR.RUNFLAG && ~ADIGATOR.EMPTYFLAG
      dbflag = cadaCheckForDerivs(b);
    else
      dbflag = 0;
    end
    bname = b.name;
    bid = b.id;
    b = parse4asgn(b);
    if strcmp(snumeric(1).type,'()') && isnumeric(x.val) && isempty(x.val)
      if nsubs > 1
        x.val = struct(snumeric(2).subs,[]);
      else
        fnames = fieldnames(b);
        structinit = cell(2,length(fnames));
        structinit(1,:) = fnames;
        x.val = struct(structinit{:});
      end
    end
    xnew.val  = subsasgn(x.val,snumeric,b);
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,asgnstring,' = ',bname,';\n']);
    end
    ADIGATOR.VARINFO.LASTOCC(bid,1) = ADIGATOR.VARINFO.COUNT;
    if PFLAG && ~ADIGATOR.EMPTYFLAG && any(ADIGATOR.VARINFO.NAMELOCS(xid,2) ~= ...
        ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.LASTOCC(:,1)==xid,2))
      if isempty(bid)
        bid = find(ADIGATOR.VARINFO.LASTOCC(:,1)==xid,1,'first');
      end
      adigatorPrintStructAsgn(b,asgnstring,bname,xid,bid);
    end
  elseif ischar(b)
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,asgnstring,' = ''',b,''';\n']);
    end
    dbflag = 0;
    xnew.val  = subsasgn(x.val,snumeric,b);
  elseif isnumeric(b) && isempty(b)
    % Empty assignment - maybe removing elements, who knows
    dbflag = 0;
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,asgnstring,' = [];\n']);
    end
    xnew.val  = subsasgn(x.val,snumeric,[]);
  end
  
  xnew.arrayflag = 1;
  x = cadastruct(xnew);
  if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 1 && ...
      (cadaflag || refasgnflag) && ~ADIGATOR.EMPTYFLAG
    if isstruct(yprior) || iscell(yprior)
      yprior = cadastruct(yprior,asgnstring,[],0);
    end
    AsgnForAsgnOvermap(x,b,yprior,cadaflag,refasgnflag);
    if cadaflag
      AsgnForAsgnInds(snumeric(1:nsubs));
    end
  end
  

  if ADIGATOR.RUNFLAG && ~ADIGATOR.EMPTYFLAG && isinf(ADIGATOR.VARINFO.NAMELOCS(ADIGATOR.VARINFO.COUNT,3)) && dbflag
    error(['variable ''',x.name,...
      ''' is either an auxiliary or global variable which was ',...
      're-assigned to have derivative information - this is not allowed.']);
  end
  ADIGATOR.VARINFO.LASTOCC([xid,xoldid],1) = xid;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;

end

end

function x = parse4asgn(x)
global ADIGATOR
if isa(x,'cadastruct') && ~x.arrayflag
 x = x.val;
end
if isa(x,'cada')
  ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
  x.id = ADIGATOR.VARINFO.COUNT;
elseif isa(x,'cadastruct')
  ADIGATOR.VARINFO.LASTOCC(x.id,1) = ADIGATOR.VARINFO.COUNT;
  x = x.val;
elseif isstruct(x) && ~isempty(x)
  fnames = fieldnames(x);
  for I = 1:length(fnames)
    F = fnames{I};
    x.(F) = parse4asgn(x.(F));
  end
end

end

function [refstr,subs,cadaflag] = parseIndex(subs)
global ADIGATOR
cadaflag = false;
if isempty(subs)
  refstr = [];
  return
end
s1 = subs{1};
if isa(s1,'cada')
  cadaflag = true;
  s1func = s1.func;
  refstr = s1func.name;
  ADIGATOR.VARINFO.LASTOCC(s1.id,1) = ADIGATOR.VARINFO.COUNT;
  if ~isempty(s1func.value)
    subs{1} = s1.func.value;
  elseif ADIGATOR.EMPTYFLAG
    subs{1} = 1;
  elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
    
  else
    error('Cannot reference with unknown index')
  end
elseif isnumeric(s1)
  if length(s1) > 1
    if ADIGATOR.PRINT.FLAG
      refstr = cadaindprint(s1);
    else
      refstr = 'holder';
    end
  elseif islogical(s1)
    if s1
      refstr = 'true';
    else
      refstr = 'false';
    end
  else
    refstr = sprintf('%1.0f',s1);
  end
elseif strcmp(s1,':')
  refstr = ':';
else
  error('?? unknown reference index')
end
if length(subs) == 2
  s2 = subs{2};
  if isa(s2,'cada')
    cadaflag = true;
    s2func = s2.func;
    refstr2 = s2func.name;
    ADIGATOR.VARINFO.LASTOCC(s2.id,1) = ADIGATOR.VARINFO.COUNT;
    if ~isempty(s2func.value)
      subs{2} = s2.func.value;
    elseif ADIGATOR.EMPTYFLAG
      subs{2} = 1;
    elseif ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 2
      
    else
      error('Cannot reference with unknown index')
    end
  elseif isnumeric(s2)
    if length(s2) > 1
      if ADIGATOR.PRINT.FLAG
        refstr2 = cadaindprint(s2);
      else
        refstr2 = 'holder';
      end
    elseif islogical(s2)
      if s2
        refstr2 = 'true';
      else
        refstr2 = 'false';
      end
    else
      refstr2 = sprintf('%1.0f',s2);
    end
  elseif strcmp(s2,':')
    refstr2 = ':';
  else
    error('?? unknown reference index')
  end
  refstr = [refstr,',',refstr2];
elseif length(subs) > 2
  error('Can only use 2D cell/structure arrays')
end

end

function count = structconvert(y,yname,count)
% need to convert a scalar structure into a structure array - need to
% modify VARINFO so that we can catch the scalar structure earlier and
% modify it.
global ADIGATOR
% "bad" cases we need to check for
% 1. Multiple assignments without () index
% v.f = f(x) 
% v.dx = f2(x)
% v(2).f = 3
% This will mess up the global count if we change the first occurance of v
%
% 2. Use of the variable prior to this catch
% v = struct('f',x); (or v.f = sin(x))
% w = v.f+1;
% w(2).f = 3;
% This too will mess up our varcount

y = orderfields(y);
fnames = fieldnames(y);
for I = 1:length(fnames)
  F = fnames{I};
  xi = y.(F);
  if isa(xi,'cada') || isa(xi,'cadastruct')
    xid = xi.id;
    
    if ADIGATOR.VARINFO.LASTOCC(xid,1) > xid
      error(['Error with variable:',yname,...
        ' -- something like: v.f = calc; w = calc(v.f); v(2).f = calc; ',...
        'please use v(1).f = calc, or instantiate to structure array']);
    elseif count > 0 && xid ~= count+1
      error(['Error with variable:',yname,...
        ' -- something like: v.f = calc; v.f2 = calc; v(2).f = calc; ',...
        'please use v(1).f = calc, v(1).f2 = calc or instantiate to structure array']);
    end
    count = xid;
    adigatorAssignImpVarNames(xid,yname,1);
  elseif isstruct(xi)
    count = structconvert(xi,yname,count);
  end
end

end

function IncreaseForAsgnCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
ADIGATORFORDATA(INNERLOC).COUNT.STRUCTASGN =...
  ADIGATORFORDATA(INNERLOC).COUNT.STRUCTASGN + 1;
end

function AsgnForAsgnOvermap(x,b,yprior,cadaflag,cadareflag)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
ASGNCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTASGN;
if cadaflag
  % Store Overmaps
  b = cadaUnionVars(b,yprior);
  if isempty(ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS)
    ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS = cell(1,2);
    ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2} = b;
    ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{1} = ...
      ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.COUNT,1);
  elseif isa(ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2},'cada') || ...
      isa(ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2},'cadastruct')
    % Is an intermediate variable coming into this, we have no idea how it
    % will change, so just change it here.
    ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2} =...
      cadaUnionVars(b,ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2});
  else
    bOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(b.id,1);
    bOver     = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2};
    if bOverLoc1 && bOver ~= bOverLoc1
      ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2} = bOverLoc1;
    end
  end
end
if cadareflag
  SUBSASGNCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.SUBSASGN; 
  % CADA subsasgn gets called prior to this - replace the overmap used for
  % the CADA subsasgn with the overmap of B
  if isempty(ADIGATORFORDATA(INNERLOC).SUBSASGN(SUBSASGNCOUNT).VARS)
    ADIGATORFORDATA(INNERLOC).SUBSASGN(SUBSASGNCOUNT).VARS = cell(1,2);
  end
  if cadaflag
    ADIGATORFORDATA(INNERLOC).SUBSASGN(SUBSASGNCOUNT).VARS{1} = ...
      ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2};
  else
    ADIGATORFORDATA(INNERLOC).SUBSASGN(SUBSASGNCOUNT).VARS{1} = ...
      cadaUnionVars(b,ADIGATORFORDATA(INNERLOC).SUBSASGN(SUBSASGNCOUNT).VARS{1});
  end
end
% emptyflag = -1 => never empty
% emptyflag = 0  => sometimes empty
% emptyflag = 1  => always empty
% ^ we only care if cadaflag is true.

if cadareflag
  NUMvod = ADIGATOR.NVAROFDIFF;
  femptyflagc   = -1;
  demptyflagc   = -ones(1,NUMvod);
  if isnumeric(yprior) && isempty(yprior)
    femptyflagc = 1;
    demptyflagc(:) = 1;
  elseif isa(yprior,'cada')
    ypfunc = yprior.func;
    ypderiv = yprior.deriv;
    if any(ypfunc.size==0)
      femptyflagc    = 1;
      demptyflagc(:) = 1;
    else
      for Vcount = 1:NUMvod
        if isempty(ypderiv(Vcount).nzlocs)
          demptyflagc(Vcount) = 1;
        end
      end
    end
  end
  emptyflagc = [femptyflagc demptyflagc];
  
  if ~isempty(ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS)
    emptyflagp = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{1}(2:end);
    emptyflag = zeros(size(emptyflagp));
    eqlocs = emptyflagp == emptyflagc;
    emptyflag(eqlocs) = emptyflagp(eqlocs);
  else
    emptyflag = emptyflagc;
  end
  ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{1} = [cadareflag,emptyflag];
else
  ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{1} = 0;
end
end

function AsgnForAsgnInds(s)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
ASGNCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTASGN;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
% Store Indices
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES = s(:);
else
  if isempty(ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES)
    ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES = ...
      struct('type',cell(0),'subs',cell(0));
  end
  ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES...
    (1:length(s),ITERCOUNT) = s;
end
end

function [y,returnflag] = ForSubsAsgn(x,s,b)
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
ASGNCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTASGN;
NUMvod    = ADIGATOR.NVAROFDIFF;
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
NDstr     = sprintf('%1.0f',ADIGATOR.DERNUMBER);

if isempty(ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS)
  % Is no CADA reference - don't have to do anything special
  returnflag = 0;
  y = [];
  return
else
  returnflag = 1;
end
bOver = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{2};
xOver = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).VARS{1};
if ~isa(xOver,'cadastruct')
 xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
end
if ~isa(bOver,'cada') && isnumeric(bOver)
  bOver = ADIGATORVARIABLESTORAGE.OVERMAP{bOver};
end

bid = b.id;

nsubs = length(s);
cadarefflag = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{1}(1);
if cadarefflag
  nsubs = nsubs-1;
end


xname = x.name;
xid   = x.id;
yid = ADIGATOR.VARINFO.COUNT;
refstr = cell(1,nsubs);
for subcount = 1:nsubs
  i = s(subcount);
  switch i.type
    case '.'
      refstr{subcount} = ['.',i.subs];
    case '()'
      refstri = parseIndex(i.subs);
      refstr{subcount} = ['(',refstri,')'];
    case '{}'
      refstri = parseIndex(i.subs);
      refstr{subcount} = ['{',refstri,'}'];
    otherwise
      error('?? invalid reference type')
  end
end
y = xOver;
y.id = yid;

if ADIGATOR.EMPTYFLAG
  ADIGATOR.VARINFO.LASTOCC([xid,bid],1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.LASTOCC(yid,1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  return
end


ynameloc = ADIGATOR.VARINFO.NAMELOCS(yid,1);
yname    = ADIGATOR.VARINFO.NAMES{ynameloc};
bnameloc = ADIGATOR.VARINFO.NAMELOCS(bid,1);
oldbnameloc = ADIGATOR.VARINFO.NAMELOCS(bid,:);
if bnameloc > 0
  boldname = ADIGATOR.VARINFO.NAMES{bnameloc}; 
end
y.name = yname;
refstr   = [yname,cell2mat(refstr)];

structrefs     = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{2};
emptyderivchks = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{3};
sizechks       = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{4};
CountName      = ADIGATORFORDATA(INNERLOC).COUNTNAME;

inputstruct.CountName  = CountName;
inputstruct.fid        = fid;
inputstruct.indent     = indent;
inputstruct.refflag    = 0;
inputstruct.NDstr      = NDstr;
if cadarefflag
  % If there is a reference s.t. s(i).f(j) = x*2
  % Need to print out bOver = s(i).f, bOver(j) = b. 
  Scount = 1;
  % Change up bnameloc to make this an intermediate variable
  oldynameloc = ADIGATOR.VARINFO.NAMELOCS(yid,:);
  if ADIGATOR.VARINFO.NAMELOCS(yid-1,1) == 0
    tempbnameloc = [0 ADIGATOR.VARINFO.NAMELOCS(yid-1,2)+1 0];
  else
    tempbnameloc = [0 1 0];
  end
  ADIGATOR.VARINFO.NAMELOCS(bid,:) = tempbnameloc;
  ADIGATOR.VARINFO.NAMES{ynameloc} = refstr;
  
  bstr = cadafuncname(bid);
  [ystr,DPFLAG] = cadafuncname();
  
  % b = y(i) where y is a structure
  bfunc = bOver.func;    bderiv = bOver.deriv;
  bMrow = bfunc.size(1); bNcol = bfunc.size(2);
  if isinf(bMrow)
    bvec = 1;
  elseif isinf(bNcol)
    bvec = 2;
  else
    vecDim = [];
    bvec = 0;
  end
  emptyfuncheck = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{1}(2);
  if bvec
    % Try to find a vecdim...
    cadabfunc = b.func;
    if isinf(cadabfunc.size(1))
      vecDim = ['size(',cadabfunc.name,',1)'];
    elseif isinf(cadabfunc.size(2))
      vecDim = ['size(',cadabfunc.name,',2)'];
    elseif bvec == 1
      vecDim = ['size(',ystr,',1)'];
    else
      vecDim = ['size(',ystr,',2)'];
    end
  end
  bfunc.name = bstr;
  
  emptycheck = emptyderivchks;
  inputstruct.emptycheck = emptycheck;
  inputstruct.vvec       = bvec;
  inputstruct.vecDim     = vecDim;
  inputstruct.refflag    = 1;
  
  
  if DPFLAG && emptyfuncheck < 1
    for Vcount = 1:NUMvod
      if ~isempty(bderiv(Vcount).nzlocs)
        emptyderivcheck = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).FLAGS{1}(2+Vcount);
        if emptyderivcheck < 1
          dbstr = cadadername(bstr,Vcount,bid);
          dystr = cadadername(ystr,Vcount);
          bderiv(Vcount).name = dbstr;
          IndName  = ...
            ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES{Vcount,Scount,1};
          IndFlags = ...
            ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES{Vcount,Scount,3};
          nzb = size(bderiv(Vcount).nzlocs,1);
          inputstruct.dxiStr   = dystr;
          inputstruct.dyStr    = dbstr;
          inputstruct.nzd      = nzb;
          inputstruct.IndName  = IndName;
          inputstruct.IndFlags = IndFlags;
          inputstruct.emptyfieldcheck = (emptyderivcheck == 0 & logical(bvec));
          cadaloopstructderivref(inputstruct);
        end
      end
    end
  end

  sizecheck = sizechks;
  
  inputstruct.sizecheck = sizecheck;
  inputstruct.yStr      = bstr;
  inputstruct.xiStr     = ystr;
  inputstruct.vsize     = [bMrow, bNcol];
  if emptyfuncheck == 1
    if bvec == 1
      fprintf(fid,[indent,bstr,' = zeros(',vecDim,',%1.0f);\n'],bNcol);
    elseif bvec == 2
      fprintf(fid,[indent,bstr,' = zeros(%1.0f,',vecDim,');\n'],bMrow);
    else
      fprintf(fid,[indent,bstr,' = zeros(%1.0f,%1.0f);\n'],bMrow,bNcol);
    end
    sizecheck = 0;
  else
    if sizecheck
      inputstruct.IndName = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).SIZES{1};
      if ~isempty(inputstruct.IndName)
        inputstruct.IndDep = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).SIZES{3}(1);
      end
      inputstruct.emptyfieldcheck = emptyfuncheck == 0 & logical(bvec);
    end
    cadaloopstructfuncref(inputstruct)
  end
  
  % should now have cada1tf1 = y(i).f
  % want to now use ynameloc for btemp
  ADIGATOR.VARINFO.NAMELOCS(bid,:) = oldbnameloc;
  ADIGATOR.VARINFO.NAMES{ynameloc} = yname;
  ADIGATOR.VARINFO.NAMELOCS(yid,:) = tempbnameloc;

  if ~sizecheck && emptyfuncheck > -1
    bfunc.emptyalloc = 1;
  end
  tempb = cada(yid,bfunc,bderiv);
  b = subsasgn(tempb,s(end),b);

  ADIGATOR.VARINFO.NAMELOCS(yid,:) = oldynameloc;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT-1;
  inputstruct.refflag    = 0;
  ADIGATOR.VARINFO.NAMELOCS(bid,:) = tempbnameloc;
end

b = cadaPrintReMap(b,bOver,bid);

if isa(b,'cadastruct')
  bname = b.name;
  if bnameloc == 0
    bnameloc = length(ADIGATOR.VARINFO.NAMES)+1;
    ADIGATOR.VARINFO.NAMELOCS(bid,:) = [bnameloc, 0 1];
  end
end

inputstruct.emptyfieldcheck = 0;
numrefvars = length(sizechks);
for Scount = 1:numrefvars
  if ~isempty(structrefs)
    % y is a structure/cell
    s2 = structrefs{Scount};
    nsubs2 = length(s2);
    refstr2 = cell(1,nsubs2);
    for subcount = 1:nsubs2
      i = s2(subcount);
      switch i.type
        case '.'
          refstr2{subcount} = ['.',i.subs];
        case '()'
          refstri = parseIndex(i.subs);
          refstr2{subcount} = ['(',refstri,')'];
        case '{}'
          refstri = parseIndex(i.subs);
          refstr2{subcount} = ['{',refstri,'}'];
        otherwise
          error('?? invalid reference type')
      end
    end
    refstr2 = cell2mat(refstr2);
    bi = subsref(b.val,s2);
  else
    refstr2 = [];
    bi = b;
  end
  
  ADIGATOR.VARINFO.NAMES{ynameloc} = [refstr,refstr2];
  if isa(b,'cadastruct')
    ADIGATOR.VARINFO.NAMES{bnameloc} = [bname,refstr2];
  end
  bstr = cadafuncname(bid);
  [ystr,DPFLAG] = cadafuncname();
  
  bfunc = bi.func; bderiv = bi.deriv;
  bMrow = bfunc.size(1); bNcol = bfunc.size(2);
  if isinf(bMrow)
    bvec = 1;
    vecDim = ['size(',bstr,',1)'];
  elseif isinf(bNcol)
    bvec = 2;
    vecDim = ['size(',bstr,',2)'];
  else
    vecDim = [];
    bvec = 0;
  end
  
  emptycheck = emptyderivchks(Scount);
  inputstruct.emptycheck = emptycheck;
  inputstruct.vvec       = bvec;
  inputstruct.vecDim     = vecDim;
  inputstruct.emptyfieldcheck = 0;
  inputstruct.Scount = Scount;
  
  if DPFLAG
    for Vcount = 1:NUMvod
      if ~isempty(bderiv(Vcount).nzlocs)
        dbstr = cadadername(bstr,Vcount,bid);
        dystr = cadadername(ystr,Vcount);
        bderiv(Vcount).name = dbstr;
        IndName  = ...
          ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES{Vcount,Scount,1};
        IndFlags = ...
          ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).INDICES{Vcount,Scount,3};
        nzb = size(bderiv(Vcount).nzlocs,1);
        inputstruct.dxiStr   = dystr;
        inputstruct.dyStr    = dbstr;
        inputstruct.nzd      = nzb;
        inputstruct.IndName  = IndName;
        inputstruct.IndFlags = IndFlags;
        cadaloopstructderivref(inputstruct);
      end
    end
  end
  
  sizecheck = sizechks(Scount);
  
  inputstruct.sizecheck = sizecheck;
  inputstruct.yStr      = bstr;
  inputstruct.xiStr     = ystr;
  inputstruct.vsize     = [bMrow, bNcol];
  if sizecheck
    inputstruct.IndName = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).SIZES{1};
    if ~isempty(inputstruct.IndName)
      inputstruct.IndDep = ADIGATORFORDATA(INNERLOC).STRUCTASGN(ASGNCOUNT).SIZES{3}(1);
    end
  end
  cadaloopstructfuncref(inputstruct)
end

ADIGATOR.VARINFO.NAMELOCS(bid,:) = oldbnameloc;
if oldbnameloc(1) > 0
  ADIGATOR.VARINFO.NAMES{oldbnameloc(1)} = boldname;
end
ADIGATOR.VARINFO.NAMES{ynameloc} = yname;
ADIGATOR.VARINFO.LASTOCC([xid,yid,bid],1) = yid;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;

end