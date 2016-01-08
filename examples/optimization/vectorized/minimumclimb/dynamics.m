function daeout = dynamics(x,u,probinfo)

CONSTANTS = probinfo.CONSTANTS;

CoF  = CONSTANTS.CoF;

h   = x(:,1);
v   = x(:,2);
fpa = x(:,3);

c1  = 392.4;    c2  = 16818;      c3  = 86.138;
c4  = 288.14;   c5  = 6.49;       c6  = 4.0519e9;
c7  = 288.08;   c8  = 5.256;      c9  = 216.64;
c10 = 9.06e8;   c11 = 1.73;       c12 = 0.157;
c13 = 6e-5;     c14 = 4.00936;    c15 = 2.2;

Nzeros = zeros(size(h));
T = Nzeros;
p = Nzeros;

ihlow     = h < 11000;
T(ihlow)  = c4 - c5*h(ihlow);
p(ihlow)  = c6*(T(ihlow)./c7).^c8;

ihhigh    = ~ihlow;
T(ihhigh) = c9;
p(ihhigh) = c10* exp(c11 - c12*h(ihhigh));

rho = c3*p./T;
q = 0.5.*rho.*v.*v.*c13;

a = c14.*sqrt(T);       M = v./a;   
Mp = cell(1,6);
for i = 1:6
  Mp{i} = M.^(i-1);
end

numeratorCD0 = Nzeros; denominatorCD0 = Nzeros; 
numeratorK   = Nzeros; denominatorK   = Nzeros;
for i = 1:6
  Mpi = Mp{i};
  if i < 6
    numeratorCD0   = numeratorCD0   + CoF(1,i).*Mpi;
    denominatorCD0 = denominatorCD0 + CoF(2,i).*Mpi;
    numeratorK     = numeratorK     + CoF(3,i).*Mpi;
  end
  denominatorK = denominatorK  + CoF(4,i).*Mpi;
end
Cd0 = numeratorCD0./denominatorCD0;
K   = numeratorK./denominatorK;
FD  = q.*(Cd0+K.*((c2^2).*(c1^2)./(q.^2)).*(u.^2));

FT = Nzeros;
for i = 1:6
  ei = Nzeros;
  for j = 1:6
    ei = ei + CoF(4+j,i).*Mp{j};
  end
  FT = FT + ei.*h.^(i-1);
end
FT = FT.*c1/c15;

hdot   = v.*sin(fpa);
vdot   = (FT-FD)./c2 - c1.*sin(fpa); 
fpadot = c1.*(u-cos(fpa))./v;   

daeout = [hdot vdot fpadot];