function y = subsref(x,s)
% Structure/cell references - we print/create a new variable for ()
% and {} array references but not for '.' field references.
% y = x(s)
global ADIGATOR
y     = x.val;
nsubs = length(s);
refstr = cell(1,nsubs);
if ~x.arrayflag
  % --------------------- Scalar Structure References ------------------- %
  for subcount = 1:nsubs
    i = s(subcount);
    if strcmp(i.type,'.')
      try
      y = y.(i.subs);
      catch
        if ADIGATOR.EMPTYFLAG
          y = [];
        end
      end
      %y = subsref(y,i);
      refstr{subcount} = i.subs;
    elseif strcmp(i.type,'()')
      y = y(1);
      refstr{subcount} = '';
    else
      dotcount = subcount - 1;
      break
    end
    if subcount == nsubs
      dotcount = nsubs;
      break
    end
    if isa(y,'cada') || isa(y,'cadastruct')
      dotcount = subcount;
      break
    end
  end
  
  if isstruct(y)
    y = cadastruct(y,[x.name,cell2mat(refstr(1:dotcount))],[],0);
  end
  
  if dotcount < nsubs
    % Call subsref on whatever is left
    y = subsref(y,s(dotcount+1:end));
  end
else
  % ------------------- Structure/Cell Array References ----------------- %
  PFLAG   = ADIGATOR.PRINT.FLAG;
  NUMvod  = ADIGATOR.NVAROFDIFF;
  fid     = ADIGATOR.PRINT.FID;
  indent  = ADIGATOR.PRINT.INDENT;
  if ADIGATOR.FORINFO.FLAG
    IncreaseForRefCount();
    if ADIGATOR.RUNFLAG == 2
      [y,returnflag] = ForSubsRef(x,s);
      if returnflag
        return
      end
      y = x.val;
    end
  end
  cadaflag = 0;
  yname = x.name;
  xid   = x.id;
  refstr = cell(1,nsubs);
  snumeric = s;
  badrefflag = 0;


  for subcount = 1:nsubs
    i = s(subcount);
    switch i.type
      case '.'
        refstr{subcount} = ['.',i.subs];
        try
        y = subsref(y,i);
        catch
          if ADIGATOR.EMPTYFLAG
            y = [];
            badrefflag = 1;
            break
          else
            refstr = [x.name,cell2mat(refstr(1:subcount))];
            error(['??? cannot find field:  ', refstr]);
          end
        end
      case '()'
        [refstri,i.subs,cadaflag2] = parseIndex(i.subs);
        if ADIGATOR.EMPTYFLAG
          try
            y = subsref(y,i);
          catch
            if ~isempty(y)
              y = y(1);
            elseif isnumeric(y) && subcount == nsubs
              % Assume this is a cada object
              subcount = subcount-1; %#ok<FXSET>
              break
            end
          end
        else
          y = subsref(y,i);
        end
        cadaflag = or(cadaflag,cadaflag2);
        refstr{subcount} = ['(',refstri,')'];
        snumeric(subcount) = i;
      case '{}'
        [refstri,i.subs,cadaflag2] = parseIndex(i.subs);
        cadaflag = or(cadaflag,cadaflag2);
        snumeric(subcount) = i;
        refstr{subcount} = ['{',refstri,'}'];
        if ADIGATOR.EMPTYFLAG
          try
            y = subsref(y,i);
          catch
            y = y{1};
          end
        else
          y = subsref(y,i);
        end
      otherwise
        error('?? invalid reference type')
    end
    if isa(y,'cada')
      break
    end
  end
  
  if badrefflag
    if strcmp(s(nsubs).type,'()')
      if ADIGATOR.DERNUMBER > 1
        subcount = nsubs-1;
      else
        error('Cannot determine structure reference due to inconsistent structure/cell array entries');
      end
    else
      subcount = nsubs;
    end
  end
  
  if isnumeric(y) && isempty(y)
    yfunc  = struct('name',[],'size',[0 0],'zerolocs',[],'value',[]);
    yderiv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
    y = cada(ADIGATOR.VARINFO.COUNT,yfunc,yderiv);
  end
  
  if PFLAG
    if subcount < nsubs
      refstr = [yname,cell2mat(refstr(1:subcount))];
    else
      refstr = [yname,cell2mat(refstr)];
    end
  end
  
  if isa(y,'cada') && subcount == nsubs
    % output is a cada variable - we are making this a new variable
    %y.id = ADIGATOR.VARINFO.COUNT;
    [funcstr,DPFLAG] = cadafuncname();
    yfunc   = y.func;
    yderiv  = y.deriv;
    yfunc.name = funcstr;
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      nameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
      oldname = ADIGATOR.VARINFO.NAMES{nameloc};
      ADIGATOR.VARINFO.NAMES{nameloc} = refstr;
      funrefstr = cadafuncname(xid);
    end

    %y.func.name = funcstr;
    for Vcount = 1:NUMvod
      if ~isempty(yderiv(Vcount).nzlocs)
        derivstr = cadadername(funcstr,Vcount);
        yderiv(Vcount).name = derivstr;
        if DPFLAG && ~ADIGATOR.EMPTYFLAG
          % Need to figure out the derivative reference..
          derrefstr = cadadername(funrefstr,Vcount,xid);
          fprintf(fid,[indent,derivstr,' = ',derrefstr,';\n']);
        end
      end
    end
    
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,funcstr,' = ',funrefstr,';\n']);
      ADIGATOR.VARINFO.NAMES{nameloc} = oldname;
    end
    y = cada(ADIGATOR.VARINFO.COUNT,yfunc,yderiv);
    ADIGATOR.VARINFO.LASTOCC([ADIGATOR.VARINFO.COUNT, xid],1) = ADIGATOR.VARINFO.COUNT;
    ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
  elseif isstruct(y) || iscell(y)
    % Make this a cadastruct no matter what it is
    yid = ADIGATOR.VARINFO.COUNT;
    if ADIGATOR.RUNFLAG == 2
      nameloc = ADIGATOR.VARINFO.NAMELOCS(yid,1);
      if nameloc > 0
        yname = ADIGATOR.VARINFO.NAMES{nameloc};
      else
        yname = sprintf(['cada',NDstr,'s%1.0f'],ADIGATOR.VARINFO.NAMELOCS(yid,2));
      end
      
    else
      yname = 'cadadummystruct';
    end
    if isinf(ADIGATOR.VARINFO.NAMELOCS(xid,3))
      ADIGATOR.VARINFO.NAMELOCS(yid,3) = ADIGATOR.VARINFO.NAMELOCS(xid,3);
    end
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      fprintf(fid,[indent,yname,' = ',refstr,';\n']);
    end
    if PFLAG && ~ADIGATOR.EMPTYFLAG && ...
        ADIGATOR.VARINFO.NAMELOCS(xid,2) ~= ADIGATOR.VARINFO.NAMELOCS(yid,2)
      adigatorPrintStructAsgn(y,yname,refstr,yid,xid);
    end
    y = cadastruct(y,yname,yid,1);
    ADIGATOR.VARINFO.LASTOCC([xid,yid],1) = yid;
    ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT + 1;
  elseif ~isa(y,'cada')
    % String, do nothing?
    cadaflag = 0;
  else
    y.id = ADIGATOR.VARINFO.COUNT;
  end
  
  if ADIGATOR.FORINFO.FLAG && ADIGATOR.RUNFLAG == 1 && ...
      (cadaflag || subcount < nsubs) && ~ADIGATOR.EMPTYFLAG
    AsgnForRefOvermap(x,y,cadaflag,subcount < nsubs);
    %count2 = ParseStructForLoop(snumeric(1:subcount),y,0);
    if cadaflag
      AsgnForRefInds(snumeric(1:subcount));
    end
  end
  
  if subcount < nsubs
    % This is a CADA subsref
    y.id = [];
    if PFLAG && ~ADIGATOR.EMPTYFLAG
      % When printing, we will modify the name of y so that CADA subsref
      % prints this out properly.
      nameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
      oldname = ADIGATOR.VARINFO.NAMES{nameloc};
      ADIGATOR.VARINFO.NAMES{nameloc} = refstr;
      % Rename y
      if isa(y,'cada')
        funcstr = cadafuncname(xid);
        yderiv = y.deriv;
        for Vcount = 1:NUMvod
          if ~isempty(yderiv(Vcount).nzlocs)
            yderiv(Vcount).name = cadadername(funcstr,Vcount,xid);
          end
        end
        y.deriv = yderiv;
        y.func.name = funcstr;
      end
      % Call subsref
      y = subsref(y,s(subcount+1:nsubs));
      % Replace VARINFO.NAMES
      ADIGATOR.VARINFO.NAMES{nameloc} = oldname;
    else
      y = subsref(y,s(subcount+1:nsubs));
    end
    ADIGATOR.VARINFO.LASTOCC(xid,1) = ADIGATOR.VARINFO.COUNT-1;
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

function IncreaseForRefCount()
global ADIGATOR ADIGATORFORDATA
INNERLOC = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTREF + 1;
ADIGATORFORDATA(INNERLOC).COUNT.STRUCTREF = REFCOUNT;
% if ADIGATOR.RUNFLAG == 1 && length(ADIGATORFORDATA(INNERLOC).STRUCTREF) ...
%     < REFCOUNT
%   if REFCOUNT == 1
%     ADIGATORFORDATA(INNERLOC).STRUCTREF = ...
%       struct('VARS',[],'FLAGS',[],'SIZES',[],'INDICES',[],'LOCS',[]);
%   end
%   ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT,1).VARS = [];
% end
end

function AsgnForRefOvermap(x,y,cadaflag,cadareflag)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTREF;
if cadaflag
% Store Overmaps
if isempty(ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS)
  ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS = cell(1,2);
  if isempty(x.id)
    ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1} = [];
  else
    OUTERLOC   = ADIGATOR.FORINFO.OUTERLOC;
    StartCount = ADIGATORFORDATA(OUTERLOC).START;
    EndCount   = ADIGATORFORDATA(OUTERLOC).END;
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    xOverLoc2 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,2);
    if xOverLoc1 && x.id >= StartCount && x.id <= EndCount
      ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1} = xOverLoc1;
    elseif xOverLoc2 && any(ADIGATOR.VARINFO.OVERMAP.FOR(StartCount:EndCount,1)==xOverLoc2)
      ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1} = xOverLoc2;
    else
      ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1} = x;
    end
  end
  ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{2}  = y;
else
  ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{2}  =...
    cadaUnionVars(ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{2},y);
  
  xOver = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1};
  if isnumeric(xOver) && ~isempty(x.id)
    xOverLoc1 = ADIGATOR.VARINFO.OVERMAP.FOR(x.id,1);
    if xOverLoc1 && xOver ~= xOverLoc1
      ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1} = xOverLoc1;
    end
  end
end
end
if cadareflag
  SUBSREFCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.SUBSREF+1;
  if isempty(ADIGATORFORDATA(INNERLOC).SUBSREF(SUBSREFCOUNT).VARS)
    ADIGATORFORDATA(INNERLOC).SUBSREF(SUBSREFCOUNT).VARS = cell(1,2);
  end
  if cadaflag
    ADIGATORFORDATA(INNERLOC).SUBSREF(SUBSREFCOUNT).VARS{1} = ...
      ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{2};
  else
    ADIGATORFORDATA(INNERLOC).SUBSREF(SUBSREFCOUNT).VARS{1} = ...
      cadaUnionVars(y,ADIGATORFORDATA(INNERLOC).SUBSREF(SUBSREFCOUNT).VARS{1});
  end
end

ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).FLAGS{1} = cadareflag;
end

function count = ParseStructForLoop(s,y,count)
if isa(y,'cadastruct')
  y = y.val;
end

if isa(y,'cada') || (isnumeric(y) && isempty(y))
  count = count+1;
  AsgnForRefInds(s,count);
elseif isstruct(y)
  Fnames = fieldnames(y);
  s2 = s;
  if numel(y) == 1
    s2(end+1).type = '.';
    for K = 1:length(Fnames)
      F = Fnames{K};
      s2(end).subs = F;
      count = ParseStructForLoop(s2,y.(F),count);
    end
  elseif numel(y) > 1
    s2(end+1).type = '()';
    s2(end+1).type = '.';
    for K = 1:length(Fnames)
      F = Fnames{K};
      for I = 1:numel(y)
        s2(end-1).subs = I;
        s2(end).subs = F;
        count = ParseStructForLoop(s2,y(I).(F),count);
      end
    end
  end
elseif iscell(y)
  s2 = s;
  s2(end+1).type = '{}';
  for I = 1:numel(y)
    s2(end).subs = I;
    count = ParseStructForLoop(s2,y{I},count);
  end
end
end

function AsgnForRefInds(s)
global ADIGATOR ADIGATORFORDATA
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTREF;
ITERCOUNT = ADIGATORFORDATA(INNERLOC).COUNT.ITERATION;
% Store Indices
if ITERCOUNT == 1
  ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).INDICES = s(:);
else
  if isempty(ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).INDICES)
    ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).INDICES = ...
      struct('type',cell(0),'subs',cell(0));
  end
  ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).INDICES...
    (1:length(s),ITERCOUNT) = s;
end
end

function [y,returnflag] = ForSubsRef(x,s)
global ADIGATOR ADIGATORFORDATA ADIGATORVARIABLESTORAGE
INNERLOC  = ADIGATOR.FORINFO.INNERLOC;
REFCOUNT  = ADIGATORFORDATA(INNERLOC).COUNT.STRUCTREF;
NUMvod    = ADIGATOR.NVAROFDIFF;
fid       = ADIGATOR.PRINT.FID;
indent    = ADIGATOR.PRINT.INDENT;
NDstr     = sprintf('%1.0f',ADIGATOR.DERNUMBER);

if isempty(ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS)
  % Is no CADA reference - don't have to do anything special
  returnflag = 0;
  y = [];
  return
else
  returnflag = 1;
end
y = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{2};
xOver = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).VARS{1};
if ~isa(xOver,'cadastruct')
 xOver = ADIGATORVARIABLESTORAGE.OVERMAP{xOver};
end
x = cadaPrintReMap(x,xOver,x.id);

%yfunc = y.func;
%yMrow = yfunc.size(1); yNcol = yfunc.size(2);


nsubs = length(s);
cadarefflag = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).FLAGS{1};
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
      [refstri,i.subs] = parseIndex(i.subs);
      refstr{subcount} = ['(',refstri,')'];
    case '{}'
      [refstri,i.subs] = parseIndex(i.subs);
      refstr{subcount} = ['{',refstri,'}'];
    otherwise
      error('?? invalid reference type')
  end
end


if ADIGATOR.EMPTYFLAG
  ADIGATOR.VARINFO.LASTOCC(xid,1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.LASTOCC(yid,1) = ADIGATOR.VARINFO.COUNT;
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
  return
end

refstr = [xname,cell2mat(refstr)];
xnameloc = ADIGATOR.VARINFO.NAMELOCS(xid,1);
xoldname = ADIGATOR.VARINFO.NAMES{xnameloc};
ynameloc = ADIGATOR.VARINFO.NAMELOCS(yid,1);
if ynameloc > 0
  yoldname = ADIGATOR.VARINFO.NAMES{ynameloc};
end

structrefs     = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).FLAGS{2};
emptyderivchks = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).FLAGS{3};
sizechks       = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).FLAGS{4};

CountName = ADIGATORFORDATA(INNERLOC).COUNTNAME;

inputstruct.CountName  = CountName;
inputstruct.fid        = fid;
inputstruct.indent     = indent;
inputstruct.refflag    = 1;
inputstruct.NDstr      = NDstr;

% Change names of y accordingly.
if isa(y,'cadastruct')
  % y is structure/cell
  if ynameloc > 0
    yname = ADIGATOR.VARINFO.NAMES{ynameloc};
  else
    yname = sprintf(['cada',NDstr,'s%1.0f'],ADIGATOR.VARINFO.NAMELOCS(yid,2));
  end
  y.id = yid;
  y.name = yname;
  if length(y.val) > 1
    TempName = 'cadats1';
    fprintf(fid,[indent,TempName,' = ',refstr,';\n']);
    refstr = TempName;
  end
else
  % y is cada
  yfunc = y.func;
  yderiv = y.deriv;
  [funcstr,DPFLAG] = cadafuncname();
  yfunc.name = funcstr;
  for Vcount = 1:NUMvod
    if ~isempty(yderiv(Vcount).nzlocs)
      derivname = cadadername(funcstr,Vcount);
      yderiv(Vcount).name = derivname;
    end
  end
  y = cada(yid,yfunc,yderiv);
end

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
    yi = subsref(y.val,s2);
  else
    refstr2 = [];
    yi = y;
  end
  % Change NAMEs up a touch.
  ADIGATOR.VARINFO.NAMES{xnameloc} = [refstr,refstr2];
  funrefstr = cadafuncname(xid);
  if ~isempty(structrefs)
    ADIGATOR.VARINFO.NAMES{ynameloc} = [yname,refstr2];
    [funcstr,DPFLAG] = cadafuncname();
  end
  
  yfunc = yi.func; yderiv = yi.deriv;
  yMrow = yfunc.size(1); yNcol = yfunc.size(2);
  TF1 = ['cada',NDstr,'tempf1'];
  if isinf(yMrow)
    yvec = 1;
    fprintf(fid,[indent,TF1,' = ',funrefstr,';\n']);
    vecDim = ['size(',TF1,',1)'];
  elseif isinf(yNcol)
    yvec = 2;
    fprintf(fid,[indent,TF1,' = ',funrefstr,';\n']);
    vecDim = ['size(',TF1,',2)'];
  else
    vecDim = [];
    yvec = 0;
  end
  
  emptycheck = emptyderivchks(Scount);
  inputstruct.emptycheck      = emptycheck;
  inputstruct.vvec            = yvec;
  inputstruct.vecDim          = vecDim;
  inputstruct.emptyfieldcheck = 0;
  inputstruct.Scount = Scount;
  
  if DPFLAG
    for Vcount = 1:NUMvod
      if ~isempty(yderiv(Vcount).nzlocs)
        derivstr = cadadername(funcstr,Vcount);
        derivrefstr = cadadername(funrefstr,Vcount,xid);
        IndName  = ...
          ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).INDICES{Vcount,Scount,1};
        IndFlags = ...
          ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).INDICES{Vcount,Scount,3};
        nzy = size(yderiv(Vcount).nzlocs,1);
        inputstruct.dxiStr   = derivrefstr;
        inputstruct.dyStr    = derivstr;
        inputstruct.nzd      = nzy;
        inputstruct.IndName  = IndName;
        inputstruct.IndFlags = IndFlags;
        cadaloopstructderivref(inputstruct);
      end
    end
  end
  
  sizecheck = sizechks(Scount);
  
  inputstruct.sizecheck = sizecheck;
  inputstruct.yStr      = funcstr;
  if logical(yvec)
    inputstruct.xiStr   = TF1;
  else
    inputstruct.xiStr   = funrefstr;
  end
  inputstruct.vsize     = [yMrow, yNcol];
  if sizecheck
    inputstruct.IndName = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).SIZES{1};
    if ~isempty(inputstruct.IndName)
      inputstruct.IndDep = ADIGATORFORDATA(INNERLOC).STRUCTREF(REFCOUNT).SIZES{3}(1);
    end
  end

  cadaloopstructfuncref(inputstruct)
end

ADIGATOR.VARINFO.NAMES{xnameloc} = xoldname;
if ynameloc > 0
  ADIGATOR.VARINFO.NAMES{ynameloc} = yoldname;
end

if cadarefflag
  y = subsref(y,s(nsubs+1));
  ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT - 1;
end

ADIGATOR.VARINFO.LASTOCC([xid,yid],1) = yid;
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
end
