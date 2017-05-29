% This file compares computing Jacobian-vector products and Hessian-vector
% products in two different ways. Jacobian-vector products may be computed
% by computing the Jacobian and then multiplying by the vector, or directly
% computing the Jacobian via seeding. Hessian-vector products may be
% computed in a similar manner.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
clc
fprintf ('AdiGator example: %s\n', mfilename ('fullpath')) ;
% Will use f(x) = arrow(x) and g(x) = arrowsum(x) = sum(f(x));
% Consider x = v*t, t scalar
n = 2^10;
x = rand(n,1);
v = rand(n,1);
opts.overwrite = 1;
opts.echo = 0;

% ------------------------- 1st Deriv of f wrt x ------------------------ %
% Want p = Jf(x)*v - can make function which computes Jf(x)
tic
x_x = adigatorCreateDerivInput([n 1],'x');
adigator('arrow',{x_x},'arrow_x',opts);
jac.gen = toc;

x1.f = x;
x1.dx = ones(n,1);
tic
for i = 1:10
  f1 = arrow_x(x1);
  J1 = sparse(f1.dx_location(:,1),f1.dx_location(:,2),f1.dx,n,n);
  p1 = J1*v;
end
jac.eval = toc/10;


% -------------------------- 1st Deriv of f wrt t ----------------------- %
% df/dt = Jf(x)*v = p
tic
derinfo.vodname = 't'; % Take derivative wrt variable called 't'
derinfo.vodsize = [1 1]; % 't' is a scalar
derinfo.nzlocs  = [(1:n).',ones(n,1)]; % input 'x' has full column deriv wrt 't'
x_t = adigatorCreateDerivInput([n 1],derinfo);
adigator('arrow',{x_t},'arrow_t',opts);
jacv.gen = toc;

x2.f = x;
x2.dt = v;
tic
for i = 1:10
  f2 = arrow_t(x2);
  p2 = f2.dt;
end
jacv.eval = toc/10;

err1 = max(abs(p1-p2)./abs(p1));
display(['n = ',num2str(n)]);
display('Times for Computing Jf(x) then Multiplying')
disp(jac);
display('Times for Computing Jf(x)*v (Jacobian vector product)')
disp(jacv);
display('Max Pct Diff of J*v both ways')
disp(err1);

% -------------------------- 2nd Deriv of g wrt x ----------------------- %
adigator('arrowsum',{x_x},'arrowsum_x',opts);
tic
ax1.f = x_x;
ax1.dx = x1.dx;
adigator('arrowsum_x',{ax1},'arrowsum_xx',opts);
hes.gen = toc;

tic
for i = 1:10
  g1 = arrowsum_xx(x1);
  H = sparse(g1.dxdx_location(:,1),g1.dxdx_location(:,2),g1.dxdx,n,n);
  q1 = H*v;
end
hes.eval = toc/10;

% ------------------------- 2nd Deriv of g wrt t ------------------------ %
% dg/dx = Gg(x); (arrowsum_x)
% d[dg/dx]/dt = Hg(x)*v;
tic
ax2.f = x_t;
ax2.dx = ones(n,1);
adigator('arrowsum_x',{ax2},'arrowsum_xt',opts);
hesv.gen = toc;

x2.dx = ones(n,1);

tic
for i = 1:10
  g2 = arrowsum_xt(x2);
  q2 = g2.dxdt;
end
hesv.eval = toc/10;

err2 = max(abs(q1-q2)./abs(q1));
display('Times for Computing Hf(x) then Multiplying')
disp(hes)
display('Times for Computing Hf(x)*v (Hessian Vector Product)')
disp(hesv)
display('Max Pct Diff of H*v both ways')
disp(err2)
