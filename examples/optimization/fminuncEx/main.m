% This function solves the brown optimisation problem using fminunc, where
% first and second derivative files are generated via
% adigatorGenFiles4Fminunc
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
clc
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
n = 2^9;
xstart = rand(n,1);
solveflag = exist('fminunc','file');


% ----------------------- Solve with 1st Derivatives -------------------- %
% setup structure for adigatorGenFiles4Fminunc
setup.numvar = n;
setup.objective = 'brownf';
setup.order = 1;
% Call adigatorGenFiles4Fminunc
tic
adifuncs = adigatorGenFiles4Fminunc(setup);
gentime = toc;

if n < 2^9
  % Dont try to solve with first derivatives if u increase n too much.
if solveflag
  tic
  options = optimoptions(@fminunc,'GradObj','on',...
    'Algorithm','trust-region','Display','iter');
  [x,fval,exitflag,output] = fminunc(adifuncs.gradient,xstart,options);
  solvetime = toc;
  disp(['1st deriv gen time + solve time: ',num2str(gentime+solvetime)]);
end
end

% ----------------------- Solve with 2nd Derivatives -------------------- %
% setup structure for adigatorGenFiles4Fminunc
setup.numvar = n;
setup.objective = 'brownf';
setup.order = 2;
% Call adigatorGenFiles4Fminunc
tic
adifuncs = adigatorGenFiles4Fminunc(setup);
gentime = toc;

if solveflag
  tic
  options = optimoptions(@fminunc,'GradObj','on','Hessian','on',...
    'Algorithm','trust-region','Display','iter');
  [x,fval,exitflag,output] = fminunc(adifuncs.hessian,xstart,options);
  solvetime = toc;
  disp(['2nd deriv gen time + solve time: ',num2str(gentime+solvetime)]);
end
