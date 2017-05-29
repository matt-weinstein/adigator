% This function solves the minimum time to climb problem by supplying first
% and second derivatives which take advantage of the vectorized nature of
% the problem. A derivative file for the dynamics is created as
% dynamics_yvect and these derivatives are used to build the constraint
% jacobian within the file conswrap.m. Also, a second derivative file for
% the dynamics is created as dynamics_yyvect, these derivatives, along with
% the first derivatives, are used to build the lagrangian hessian within
% the file laghesswrap.m.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
solveflag = exist('fmincon','file');
if solveflag
options = optimset('Algorithm','interior-point','MaxFunEvals',50000,...
  'GradObj','on','GradConstr','on','Display','iter','MaxIter',1000);
end
time = zeros(4,1);
numintervals = [10,20,40];

[probinfo,upperbound,lowerbound,guess,tau] = setupproblem(numintervals(1));

% ------------------- Create Dynamics Derivative Files ------------------ %
% Since we are doing this vectorized, we only need to take derivatives once
% Inputs to dynamics file are X,U,probinfo
% X is size N by n
% U is size N by m
% Let y(t) = [x(t) u(t)], Y = [X U], taking derivatives of dynamics file
% wrt y(t).
% Note: our decision vector, z, is of form z = [Y(:);tf], and dynamics file
% is not a function of tf.
tic
n = 3; m = 1; 
% First Derivative
gX = adigatorCreateDerivInput([Inf n],...
  struct('vodname','Y','vodsize',[Inf (m+n)],...
  'nzlocs',[(1:n).' (1:n).']));
gU = adigatorCreateDerivInput([Inf m],...
  struct('vodname','Y','vodsize',[Inf (m+n)],...
  'nzlocs',[(1:m).' (n+1:m+n).']));
gOut = adigator('dynamics',{gX,gU,probinfo},'dynamics_yvect',...
  adigatorOptions('overwrite',1));
% We can now extract sparsity pattern of Fy(t)
iFy = gOut{1}.deriv.nzlocs(:,1);
jFy = gOut{1}.deriv.nzlocs(:,2);

% Second Derivative
gX  = struct('f',gX,'dY',adigatorCreateAuxInput([Inf n]));
gU  = struct('f',gU,'dY',adigatorCreateAuxInput([Inf m]));
gOut2 = adigator('dynamics_yvect',{gX,gU,probinfo},'dynamics_yyvect',...
  adigatorOptions('overwrite',1,'comments',0));
% We can now extract sparsity pattern of Fyy(t)
iFdyy = gOut2{1}.dY.deriv.nzlocs(:,1); % Note, these are the deriv locations
kFyy  = gOut2{1}.dY.deriv.nzlocs(:,2); % of F.dy which is a vector of non-zeros
% We need to switch so that we have 3 indices since Fyy(t) is of 
% dim n by (n+m) by (n+m). That is, the first index, iFdyy, is the function
% location of F.dy and kFyy is the variable location (our third location).
% We can create our first two by referencing off of iFy and jFy
% respectively.
iFyy = iFy(iFdyy);
jFyy = jFy(iFdyy);

time(1) = toc;
nummesh = 3;
for i = 1:nummesh
tic;
% ---------------------- Set Up the Problem ----------------------------- %
if i > 1
% If on second or third mesh, use previous solution as initial guess
[probinfo,upperbound,lowerbound,guess,tau] = ...
  setupproblem(numintervals(i),X,U,fval,tau);
end
probinfo.vectders = 1;
N = probinfo.N;


% ----------------------- Project Deriv Indices ------------------------- %
% We have already created the derivative file, but now that we have fixed
% the vectorized dimension, we can find the unrolled locations of the
% unrolled derivative. We do this here rather than on each derivative
% evaluation to save time.
% 1st Deriv
[probinfo.iFy, probinfo.jFy] = adigatorProjectVectLocs(N,iFy,jFy);
% 2nd Deriv
[probinfo.iFyy, JFyy,KFyy] = adigatorProjectVectLocs(N,iFyy,jFyy,kFyy);
% We want to unroll our second two dimensions
probinfo.jFyy = sub2ind([N*(n+m) N*(n+m)],JFyy,KFyy);

if solveflag
options = optimset(options,'Hessian','user-supplied','HessFcn',...
  @(z,lambda)laghesswrap(z,lambda,probinfo));
end

% --------------------------- Call fmincon ------------------------------ %
if solveflag
[z,fval] = ...
  fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
  lowerbound,upperbound,@(x)conswrap(x,probinfo),options);
else
  % just run this for 30 jacobian/hessian evaluations
  z = guess;
  ldummy.eqnonlin = ones((N-1)*n,1);
  for k = 1:30
    H = laghesswrap(z,ldummy,probinfo);
    [C,Ceq,JC,JCeq] = conswrap(z,probinfo);
    [fval,G] = basic_objwrap(z,probinfo);
  end
end

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*z(end);
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = toc;

% ---------------------------- Plot Solution ---------------------------- %
if solveflag
disp(z(end))
figure(i);
subplot(2,1,1)
plot(t,X,'-o');
xlabel('time')
ylabel('states')
title(sprintf(['Mesh %1.0f Solution in ',num2str(time(i)),'s'],i))
subplot(2,1,2);
plot(t,U,'-o');
xlabel('time')
ylabel('control')
end

end


fprintf(['Total Time Supplying Second Derivatives (vectorized): ',...
  num2str(sum(time)),'\n']);
