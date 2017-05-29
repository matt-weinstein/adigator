% This file uses both MATLAB finite differences as well as adigator in order
% to solve the brusselator ODE using ode15s.
% 
% For this problem, in general, compressed finite differencing seems to be
% the best option as far as time is concerned. However, we compute the
% sparsity matrix using adigator. You can play with N and tspan, but be
% weary that N > 100 will make the non-compressed finite differencing
% extremely slow.
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
% ----------------------- Problem Set Up -------------------------------- %
fdflag = 1;
N = 2^3;
if N > 100
  fdflag = 0;
end
tspan = [0 100];
y0 = [1+sin((2*pi/(N+1))*(1:N)); repmat(3,1,N)];

% --------------- Generate ADiGator Derivative File --------------------- %
tic
gt = adigatorCreateAuxInput([1 1]); % aux input for t
gy = adigatorCreateDerivInput([2*N, 1],'y'); % deriv input for y
gout = adigatorGenJacFile('mybrussode',{gt,gy,N});
adigatorgentime = toc;

% Get Jacobian sparsity pattern from output of adigator
S = gout.JacobianStructure;

% ---------------- Solve ODE Using Finite Difference -------------------- %
if fdflag
  options = odeset('Vectorized','on','Stats','on','RelTol',sqrt(eps),'AbsTol',sqrt(eps)); % Set options
  tic;
  [t1,y1] = ode15s(@(t,y)mybrussode(t,y,N),tspan,y0,options);   % Solve ODE
  FDtime = toc;
end
display(' ');
% ---------------- Solve ODE Using Compressed Finite Difference --------- %
options = odeset('Vectorized','on','JPattern',S,'Stats','on','RelTol',sqrt(eps),'AbsTol',sqrt(eps)); % Set options
tic;
[t2,y2] = ode15s(@(t,y)mybrussode(t,y,N),tspan,y0,options);     % Solve ODE
CFDtime = toc;
display(' ');

% --------------------- Solve ODE Using ADiGator ------------------------ %
tic
options = odeset('Vectorized','on','Jacobian',@(t,y)mybrussode_Jac(t,y,N),'JPattern',S,'Stats','on','RelTol',sqrt(eps),'AbsTol',sqrt(eps));
[t,y] = ode15s(@(t,y)mybrussode(t,y,N),tspan,y0,options);     % Solve ODE
adigatortime = toc;
display(' ');

display(['ADiGator deriv file generation time:               ',num2str(adigatorgentime)]);
display(['ODE solve time using ADiGator:                     ',num2str(adigatortime)]);
if fdflag
  display(['ODE solve time using Finite Difference:            ',num2str(FDtime)]);
end
display(['ODE solve time using Compressed Finite Difference: ',num2str(CFDtime)]);
