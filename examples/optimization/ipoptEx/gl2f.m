% c     **********
% c
% c     Subroutine dgl2fc
% c
% c     This subroutine computes the function and gradient of the
% c     Ginzburg-Landau (2-dimensional) superconductivity problem.
% c
% c     The subroutine statement is
% c
% c       subroutine dgl2fc(nx,ny,x,y,vpotx,vpoty,f,
% c    +                  gradx,grady,gradax,graday,task,vornum)
% c
% c     where
% c
% c       nx is an integer variable.
% c         On entry nx is the number of grid points in the first
% c            coordinate direction.
% c         On exit nx is unchanged.
% c
% c       ny is an integer variable.
% c         On entry ny is the number of grid points in the second
% c            coordinate direction.
% c         On exit ny is unchanged.
% c
% c       x is a double precision array of dimension nx*ny.
% c         On entry x specifies the real part of the order parameter
% c            if task = 'F', 'G', or 'FG'.
% c            Otherwise x need not be specified.
% c         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
% c            x is set according to task.
% c
% c       y is a double precision array of dimension nx*ny.
% c         On entry y specifies the imaginary part of the order parameter
% c            if task = 'F', 'G', or 'FG'.
% c            Otherwise y need not be specified.
% c         On exit y is unchanged if task = 'F', 'G', or 'FG'. Otherwise
% c            y is set according to task.
% c
% c       vpotx is a double precision array of dimension nx*ny.
% c         On entry vpotx specifies the x component of the vector
% c            potential if task = 'F', 'G', or 'FG'.
% c            Otherwise vpotx need not be specified.
% c         On exit vpotx is unchanged if task = 'F', 'G', or 'FG'.
% c            Otherwise vpotx is set according to task.
% c
% c       vpoty is a double precision array of dimension nx*ny.
% c         On entry vpoty specifies the y component of the vector
% c            potential if task = 'F', 'G', or 'FG'.
% c            Otherwise vpoty need not be specified.
% c         On exit vpoty is unchanged if task = 'F', 'G', or 'FG'.
% c            Otherwise vpoty is set according to task.
% c
% c       f is a double precision variable.
% c         On entry f need not be specified.
% c         On exit f is set to the function evaluated at x if task = 'F'
% c            or 'FG'.
% c
% c       gradx is a double precision array of dimension nx*ny.
% c         On entry gradx need not be specified.
% c         On exit gradx contains the gradient with respect to x
% c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
% c
% c       grady is a double precision array of dimension nx*ny.
% c         On entry grady need not be specified.
% c         On exit grady contains the gradient with respect to y
% c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
% c            if task = 'G' or 'FG'.
% c
% c       gradax is a double precision array of dimension nx*ny.
% c         On entry gradax need not be specified.
% c         On exit gradax contains the gradient with respect to vpotx
% c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
% c            if task = 'G' or 'FG'.
% c
% c       graday is a double precision array of dimension nx*ny.
% c         On entry graday need not be specified.
% c         On exit graday contains the gradient with respect to vpoty
% c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
% c            if task = 'G' or 'FG'.
% c
% c       task is a character*60 variable.
% c         On entry task specifies the action of the subroutine:
% c
% c            task               action
% c            ----               ------
% c             'F'     Evaluate the function at (x,y,vpotx,vpoty).
% c             'G'     Evaluate the gradient at (x,y,vpotx,vpoty).
% c             'FG'    Evaluate the function and the gradient at
% c                         (x,y,vpotx,vpoty).
% c             'XS'    Set (x,y,vpotx,vpoty) to the standard starting
% c                         point xs.
% c
% c         On exit task is unchanged.
% c
% c       vornum is an integer variable.
% c         On entry vornum is the number of vortices.
% c         On exit vornum is unchanged.
% c
% c     MINPACK-2 Project. March 1999.
% c     Argonne National Laboratory and University of Minnesota.
% c     Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright.
% c
% c     **********
function f = gl2f(z,auxdata)

nx = auxdata.nx;
ny = auxdata.ny;
hx = auxdata.hx;
hy = auxdata.hy;
vornum = auxdata.vornum;
tkappa = auxdata.vornum;

sqn  = nx*ny;
nxp1 = nx+1;
nyp1 = ny+1;

Z = reshape(z,nx*ny,4);


znx = zeros(nx,1); znyp1 = zeros(1,nyp1);
x     = [reshape(Z(:,1),nx,ny) znx; znyp1];
y     = [reshape(Z(:,2),nx,ny) znx; znyp1];
vpotx = [reshape(Z(:,3),nx,ny) znx; znyp1];
vpoty = [reshape(Z(:,4),nx,ny) znx; znyp1];

%     Enforce vortex constraint and boundary conditions.

%     Right face for order parameter and vector potential.
arg = 2*pi*vornum*(0:ny)/ny;
x(nxp1,:) = x(1,:).*cos(arg) - y(1,:).*sin(arg);
y(nxp1,:) = x(1,:).*sin(arg) + y(1,:).*cos(arg);
vpotx(nxp1,:) = vpotx(1,:);
vpoty(nxp1,:) = vpoty(1,:) + 2*pi*vornum/(ny*hy);

%     Top face for order parameter and vector potential.
x(:,nyp1) = x(:,1);
y(:,nyp1) = y(:,1);
vpotx(:,nyp1) = vpotx(:,1);
vpoty(:,nyp1) = vpoty(:,1);

% Get Common Reference Variables
i = 1:nx;             j = 1:ny;
xij = x(i,j);         yij= y(i,j);
vpotxij = vpotx(i,j); vpotyij = vpoty(i,j);

%        Compute the Condensation Energy Density
delsq = xij.^2 + yij.^2;
delsq = delsq(:);
fcond = sum(- delsq + (delsq.^2)/2)/sqn;

%        Compute the Kinetic Energy Density.
fkin = 0;
x1 = x(i+1,j) - xij.*cos(hx*vpotxij) +...
  +             yij.*sin(hx*vpotxij);
x2 = y(i+1,j) - yij.*cos(hx*vpotxij) -...
  +             xij.*sin(hx*vpotxij);
x1 = x1(:); x2 = x2(:);
fkin = fkin + sum((x1.^2+x2.^2)/(hx.^2));
x1 = x(i,j+1) - xij.*cos(hy*vpotyij) +...
  +             yij.*sin(hy*vpotyij);
x2 = y(i,j+1) - yij.*cos(hy*vpotyij) -...
  +             xij.*sin(hy*vpotyij);
x1 = x1(:); x2 = x2(:);
fkin = fkin + sum((x1.^2+x2.^2)/(hy.^2));
fkin = fkin/sqn;

%        Compute the Magnetic Field Energy Density.
xy = (vpoty(i+1,j)-vpotyij)/hx -...
  +              (vpotx(i,j+1)-vpotxij)/hy;
xy = xy(:);
ffield = sum(xy.^2)*(tkappa.^2)/sqn;

%  Compute f
f = fcond + fkin + ffield;