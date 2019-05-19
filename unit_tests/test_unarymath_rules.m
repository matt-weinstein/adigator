function violations = test_unarymath_rules()

mkdir('tmp');
addpath('tmp');
fid = fopen('cadaunarymath.m','r');
frewind(fid);

l = fgetl(fid);
while ~isnumeric(l)
  if strcmp(strtrim(l),'function dydx = getdydx(x,callerstr)')
    break;
  end
  l = fgetl(fid);
end

fnames = cell(1,0);
while ~isnumeric(l)
  l = strtrim(l);
  if strncmp(l,'case',4)
    tmp = strtrim(l(5:end));
    fnames{end+1} = tmp(2:end-1); %#ok<AGROW>
  end
  l = fgetl(fid);
end
bad = zeros(size(fnames));
for ii = 1:length(fnames)
  bad(ii) = test_this(fnames{ii});
end

for ii = 1:length(fnames)
  if bad(ii)
    fprintf('%s - %d violations\n',fnames{ii},bad(ii));
  end
end
violations = sum(bad);

end

function bad = test_this(fun)
fname = ['test_',fun];
ffname = ['tmp/',fname];
fid2 = fopen([ffname,'.m'],'w+');
fprintf(fid2,'function y = %s(x)\n',fname);
fprintf(fid2,'y = %s(x);\n',fun);
fprintf(fid2,'end\n');
fclose(fid2);

dname = [fname,'_dx'];
rehash;
ax = adigatorCreateDerivInput([1 1],'x');

adigator(fname,{ax},dname,adigatorOptions('overwrite',1));
movefile([dname,'.*'],'tmp/');
bad = 0;

xtest = [linspace(-0.9,0.9,10) linspace(-2*pi,2*pi,10) linspace(-360,360,10)];
for ii = 1:length(xtest)
  xx.f  = xtest(ii);
  xx.dx = 1;
  yy = feval(dname,xx);
  
  ee = 1e-6;
  f1 = feval(fname,xx.f);
  f2 = feval(fname,xx.f+ee);
  df = (f2-f1)/ee;
  if abs(yy.dx-df)/(1+abs(df)) > 1e-4
    bad = bad + 1;
    if isnan(f1) || isinf(f1) || abs(df) > 1e8
      % Neglect these corner cases - finite difference as wrong as AD
      bad = bad-1;
    else
      %keyboard
    end
  end
end

end