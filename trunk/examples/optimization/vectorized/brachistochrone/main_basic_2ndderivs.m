% This function solves the brachistochrone problem by supplying first and
% second derivatives without taking advantage of the vectorized nature of
% the problem. A wrapper file for the constraint and objective functions
% (which call derivative files) are created and supplied in
% basic_conswrap.m and basic_objwrap.m, respectively. The lagrangian
% hessian is computed by taking the derivaive of the file basic_laggrad
% which builds the lagrangian gradient. A wrapper file for the lagrangian
% hessian is supplied in basic_laghess.m
% 
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;

solveflag = exist('fmincon','file');
if solveflag
  options = optimset('Algorithm','interior-point','MaxFunEvals',50000,...
    'GradObj','on','GradConstr','on','Display','iter');
end

numintervals = [5,10,20,40];
time = zeros(length(numintervals),1);
for i = 1:length(numintervals)
start = tic;
% ---------------------- Set Up the Problem ----------------------------- %
if i == 1
[probinfo,upperbound,lowerbound,guess,tau] = setupproblem(numintervals(i));
else
% If on second or third mesh, use previous solution as initial guess
[probinfo,upperbound,lowerbound,guess,tau] = ...
  setupproblem(numintervals(i),X,U,fval,tau);
end
% Need to redo this optimset each time as probinfo changes and is input to
% basic_laggradwrap
if solveflag
  options = optimset(options,'Hessian','user-supplied','HessFcn',...
    @(z,lambda)basic_laggradwrap(z,lambda,probinfo));
end

% ----------------- Create Constraint Derivative File ------------------- %
% NOTE: fmincon is given the function basic_conswrap which builds the
% Jacobian by calling basic_cons_z
gz = adigatorCreateDerivInput([length(guess), 1],'z');
outputs = adigator('basic_cons',{gz,probinfo},'basic_cons_z',adigatorOptions('overwrite',1));

% -------------- Create Lagrangian Gradient Derivative File ------------- %
% NOTE: fmincon is given the function basic_laggradwrap which builds the
% Lagrangian Hessian by calling basic_laggrad_z
glambda = adigatorCreateAuxInput(outputs{2}.func.size);
gH1 = adigator('basic_laggrad',{gz,glambda,probinfo},'basic_laggrad_z',...
  adigatorOptions('overwrite',1,'comments',0));

% --------------------------- Call fmincon ------------------------------ %
if solveflag
  [z,fval] = ...
    fmincon(@(x)basic_objwrap(x,probinfo),guess,[],[],[],[],...
    lowerbound,upperbound,@(x)basic_conswrap(x,probinfo),options);
else
  % No optimization toolbox, just test that functions work
  z = guess;
  [fval,G] = basic_objwrap(z,probinfo);
  [C,Ceq,JC,JCeq] = basic_conswrap(z,probinfo);
  lambda.eqnonlin = ones(length(Ceq),1);
  [H,G2] = basic_laggradwrap(z,lambda,probinfo);
end

% --------------------------- Extract Solution -------------------------- %
t = tau(:).*fval;
X = z(probinfo.xind);
U = z(probinfo.uind);
time(i) = toc(start);

% ---------------------------- Plot Solution ---------------------------- %
figure(i);
subplot(2,1,1)
plot(t,X,'-o');
xlabel('time')
ylabel('states')
legend('x-position','y-position','speed','location','northwest')
title(sprintf(['Mesh %1.0f Solution in ',num2str(time(i)),'s'],i))
subplot(2,1,2);
plot(t,U,'-o');
xlabel('time')
ylabel('control')
end

fprintf(['Total Time Supplying Second Derivatives (non-vectorized): ',...
  num2str(sum(time)),'\n']);
