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
% c         On entry vpotx specifies the x comp1nt of the vector
% c            potential if task = 'F', 'G', or 'FG'.
% c            Otherwise vpotx need not be specified.
% c         On exit vpotx is unchanged if task = 'F', 'G', or 'FG'.
% c            Otherwise vpotx is set according to task.
% c
% c       vpoty is a double precision array of dimension nx*ny.
% c         On entry vpoty specifies the y comp1nt of the vector
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
function z = gl2st(auxdata)
% NOTE: in this 1, x,y,vpotx, vpoty are all size nx by ny
nx = auxdata.nx;
ny = auxdata.ny;
vornum = auxdata.vornum;

tkappa = 5;
hx = sqrt(vornum/2)*3/nx;
hy = sqrt(vornum/2)*3*sqrt(3)/ny;
sqn = nx*ny;
pi = 4*atan(1);
bave = 2*pi*vornum*tkappa/(sqn*hx*hy);
sqrtv = sqrt(vornum)*pi;

% Initialize z
x     = zeros(nx,ny);
y     = zeros(nx,ny);
vpotx = zeros(nx,ny);
vpoty = zeros(nx,ny);
for j = 1:ny
  ypt = (j-1)*hy;
  for i = 1:nx
    xpt = (i-1)*hx;
    x(i,j) = 1 - (sin(sqrtv*xpt/(2*3))*...
      +      sin(sqrtv*ypt/(2*sqrt(3)*3)))^2;
    vpoty(i,j) = bave*xpt/tkappa;
  end
end

z = [x(:); y(:); vpotx(:); vpoty(:)];