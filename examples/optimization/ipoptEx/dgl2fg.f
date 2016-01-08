      subroutine dgl2fg(nx,ny,x,f,fgrad,task,w,vornum)
      character*60 task
      integer nx, ny, vornum
      double precision f
      double precision x(4*nx*ny), fgrad(4*nx*ny), w(4*(nx+1)*(ny+1))
c     **********
c
c     Subroutine dgl2fg
c
c     This subroutine computes the function and gradient of the
c     Ginzburg-Landau (2-dimensional) superconductivity problem.
c
c     The subroutine statement is
c
c       subroutine dgl2fg(nx,ny,x,f,fgrad,task,w,vornum)
c
c     where
c
c       nx is an integer variable.
c         On entry nx is the number of grid points in the first
c            coordinate direction.
c         On exit nx is unchanged.
c
c       ny is an integer variable.
c         On entry ny is the number of grid points in the second
c            coordinate direction.
c         On exit ny is unchanged.
c
c       x is a double precision array of dimension 4*nx*ny.
c         On entry x specifies the vector x if task = 'F', 'G', or 'FG'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
c            x is set according to task.
c
c       f is a double precision variable.
c         On entry f need not be specified.
c         On exit f is set to the function evaluated at x if task = 'F'
c            or 'FG'.
c
c       fgrad is a double precision array of dimension 4*nx*ny.
c         On entry fgrad need not be specified.
c         On exit fgrad contains the gradient evaluated at x if
c            task = 'G' or 'FG'.
c
c       task is a character*60 variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at x.
c             'G'     Evaluate the gradient at x.
c             'FG'    Evaluate the function and the gradient at x.
c             'XS'    Set x to the standard starting point xs.
c
c         On exit task is unchanged.
c
c       w is a double precision work array of dimension 4*(nx+1)(ny+1).
c
c       vornum is an integer variable.
c         On entry vornum specifies the number of vortices.
c         On exit vornum is unchanged.
c
c     Subprograms called
c
c       MINPACK-supplied   ...   dgl2fc
c
c     MINPACK-2 Project. August 1999.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright.
c
c     **********
      double precision zero
      parameter (zero=0.0d0)

      integer ctr, i, itemp, j, k

      external dgl2fc

      itemp = (nx+1)*(ny+1)

c     Pack work array.

      if (task(1:2) .ne. 'XS') then
         ctr = 1
         do j = 1, ny
            do i = 1, nx
               k = (j-1)*nx + i
               w(ctr) = x(k)
               w(itemp+ctr) = x(nx*ny+k)
               w(2*itemp+ctr) = x(2*nx*ny+k)
               w(3*itemp+ctr) = x(3*nx*ny+k)
               ctr = ctr + 1
            end do
            ctr = ctr + 1
         end do
      end if 

      call dgl2fc(nx,ny,w(1),w(itemp+1),w(2*itemp+1),w(3*itemp+1),f,
     +            fgrad(1),fgrad(nx*ny+1),fgrad(2*nx*ny+1),
     +            fgrad(3*nx*ny+1),task,vornum)

c     Unpack work array

      if (task(1:2) .eq. 'XS') then
         ctr = 1
         do j = 1, ny
            do i = 1, nx
               k = (j-1)*nx + i
               x(k) = w(ctr)
               x(nx*ny+k) = w(itemp+ctr)
               x(2*nx*ny+k) = w(2*itemp+ctr)
               x(3*nx*ny+k) = w(3*itemp+ctr)
               ctr = ctr + 1
            end do
            ctr = ctr + 1
         end do
      end if

      end

      subroutine dgl2fc(nx,ny,x,y,vpotx,vpoty,f,gradx,grady,gradax,
     +                 graday,task,vornum)
      character*60 task
      integer nx, ny, vornum
      double precision x(nx+1,ny+1), y(nx+1,ny+1), vpotx(nx+1,ny+1),
     +                 vpoty(nx+1,ny+1), gradx(nx,ny), grady(nx,ny),
     +                 gradax(nx,ny), graday(nx,ny)
c     **********
c
c     Subroutine dgl2fc
c
c     This subroutine computes the function and gradient of the
c     Ginzburg-Landau (2-dimensional) superconductivity problem.
c
c     The subroutine statement is
c
c       subroutine dgl2fc(nx,ny,x,y,vpotx,vpoty,f,
c    +                  gradx,grady,gradax,graday,task,vornum)
c
c     where
c
c       nx is an integer variable.
c         On entry nx is the number of grid points in the first
c            coordinate direction.
c         On exit nx is unchanged.
c
c       ny is an integer variable.
c         On entry ny is the number of grid points in the second
c            coordinate direction.
c         On exit ny is unchanged.
c
c       x is a double precision array of dimension nx*ny.
c         On entry x specifies the real part of the order parameter
c            if task = 'F', 'G', or 'FG'.
c            Otherwise x need not be specified.
c         On exit x is unchanged if task = 'F', 'G', or 'FG'. Otherwise
c            x is set according to task.
c
c       y is a double precision array of dimension nx*ny.
c         On entry y specifies the imaginary part of the order parameter
c            if task = 'F', 'G', or 'FG'.
c            Otherwise y need not be specified.
c         On exit y is unchanged if task = 'F', 'G', or 'FG'. Otherwise
c            y is set according to task.
c
c       vpotx is a double precision array of dimension nx*ny.
c         On entry vpotx specifies the x component of the vector
c            potential if task = 'F', 'G', or 'FG'.
c            Otherwise vpotx need not be specified.
c         On exit vpotx is unchanged if task = 'F', 'G', or 'FG'.
c            Otherwise vpotx is set according to task.
c
c       vpoty is a double precision array of dimension nx*ny.
c         On entry vpoty specifies the y component of the vector
c            potential if task = 'F', 'G', or 'FG'.
c            Otherwise vpoty need not be specified.
c         On exit vpoty is unchanged if task = 'F', 'G', or 'FG'.
c            Otherwise vpoty is set according to task.
c
c       f is a double precision variable.
c         On entry f need not be specified.
c         On exit f is set to the function evaluated at x if task = 'F'
c            or 'FG'.
c
c       gradx is a double precision array of dimension nx*ny.
c         On entry gradx need not be specified.
c         On exit gradx contains the gradient with respect to x
c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
c
c       grady is a double precision array of dimension nx*ny.
c         On entry grady need not be specified.
c         On exit grady contains the gradient with respect to y
c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
c            if task = 'G' or 'FG'.
c
c       gradax is a double precision array of dimension nx*ny.
c         On entry gradax need not be specified.
c         On exit gradax contains the gradient with respect to vpotx
c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
c            if task = 'G' or 'FG'.
c
c       graday is a double precision array of dimension nx*ny.
c         On entry graday need not be specified.
c         On exit graday contains the gradient with respect to vpoty
c            of f evaluated at (x,y,vpotx,vpoty) if task = 'G' or 'FG'.
c            if task = 'G' or 'FG'.
c
c       task is a character*60 variable.
c         On entry task specifies the action of the subroutine:
c
c            task               action
c            ----               ------
c             'F'     Evaluate the function at (x,y,vpotx,vpoty).
c             'G'     Evaluate the gradient at (x,y,vpotx,vpoty).
c             'FG'    Evaluate the function and the gradient at
c                         (x,y,vpotx,vpoty).
c             'XS'    Set (x,y,vpotx,vpoty) to the standard starting
c                         point xs.
c
c         On exit task is unchanged.
c
c       vornum is an integer variable.
c         On entry vornum is the number of vortices.
c         On exit vornum is unchanged.
c
c     MINPACK-2 Project. March 1999.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Paul L. Plassmann, and Stephen J. Wright.
c
c     **********
      double precision five, four, one, three, two, zero
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0,
     +          five=5.0d0)

      integer i, j
      double precision arg, bave, cfac, delsq, f, fcond, ffield, fkin,
     +                 fkinx1, fkinx2, fkiny1, fkiny2, hx, hy, pi, sfac,
     +                 sqn, sqrtv, tkappa, x1, x2, xpt, xy, ypt

c     Initialize.

      tkappa = five
      hx = sqrt(vornum/two)*three/dble(nx)
      hy = sqrt(vornum/two)*three*sqrt(three)/dble(ny)
      sqn = dble(nx*ny)
      pi = four*atan(one)
      bave = two*pi*vornum*tkappa/(sqn*hx*hy)
      sqrtv = sqrt(dble(vornum))*pi

      if (task(1:2) .eq. 'XS') then

c        Initial Order Parameter.

         do 20 j = 1, ny + 1
            ypt = (dble(j)-one)*hy
            do 10 i = 1, nx + 1
               xpt = (dble(i)-one)*hx
               x(i,j) = one - (sin(sqrtv*xpt/(two*three))*
     +                  sin(sqrtv*ypt/(two*sqrt(three)*three)))**2
               y(i,j) = zero
   10       continue
   20    continue

c        Initial Vector Potential.

         do 40 j = 1, ny + 1
            do 30 i = 1, nx + 1
               xpt = (dble(i)-one)*hx
               vpotx(i,j) = zero
               vpoty(i,j) = bave*xpt/tkappa
   30       continue
   40    continue

         return

      end if

c     Enforce vortex constraint and boundary conditions.

c     Right face for order parameter and vector potential.

      do 50 j = 1, ny + 1
         arg = two*pi*vornum*(dble(j)-one)/dble(ny)
         x(nx+1,j) = x(1,j)*cos(arg) - y(1,j)*sin(arg)
         y(nx+1,j) = x(1,j)*sin(arg) + y(1,j)*cos(arg)
         vpotx(nx+1,j) = vpotx(1,j)
         vpoty(nx+1,j) = vpoty(1,j) + two*pi*vornum/(dble(ny)*hy)
   50 continue

c     Top face for order parameter and vector potential.

      do 60 i = 1, nx + 1
         x(i,ny+1) = x(i,1)
         y(i,ny+1) = y(i,1)
         vpotx(i,ny+1) = vpotx(i,1)
         vpoty(i,ny+1) = vpoty(i,1)
   60 continue

      if (task(1:1) .eq. 'F' .or. task(1:2) .eq. 'FG') then

c        Compute the Condensation Energy Density

         fcond = zero
         do 80 i = 1, nx
            do 70 j = 1, ny
               delsq = x(i,j)**2 + y(i,j)**2
               fcond = fcond - delsq + (delsq**2)/two
   70       continue
   80    continue
         fcond = fcond/sqn

c        Compute the Kinetic Energy Density.

         fkin = zero
         do 100 i = 1, nx
            do 90 j = 1, ny
               x1 = x(i+1,j) - x(i,j)*cos(hx*vpotx(i,j)) +
     +              y(i,j)*sin(hx*vpotx(i,j))
               x2 = y(i+1,j) - y(i,j)*cos(hx*vpotx(i,j)) -
     +              x(i,j)*sin(hx*vpotx(i,j))
               fkin = fkin + (x1**2+x2**2)/(hx**2)
               x1 = x(i,j+1) - x(i,j)*cos(hy*vpoty(i,j)) +
     +              y(i,j)*sin(hy*vpoty(i,j))
               x2 = y(i,j+1) - y(i,j)*cos(hy*vpoty(i,j)) -
     +              x(i,j)*sin(hy*vpoty(i,j))
               fkin = fkin + (x1**2+x2**2)/(hy**2)
   90       continue
  100    continue
         fkin = fkin/sqn

c        Compute the Magnetic Field Energy Density.

         ffield = zero
         do 120 i = 1, nx
            do 110 j = 1, ny
               xy = (vpoty(i+1,j)-vpoty(i,j))/hx -
     +              (vpotx(i,j+1)-vpotx(i,j))/hy
               ffield = ffield + xy**2
  110       continue
  120    continue
         ffield = ffield*(tkappa**2)/sqn
         f = fcond + fkin + ffield
      end if

      if (task(1:1) .eq. 'G' .or. task(1:2) .eq. 'FG') then

         do 140 j = 1, ny
            do 130 i = 1, nx
               gradx(i,j) = x(i,j)*(-one+x(i,j)**2+y(i,j)**2)
               gradx(i,j) = gradx(i,j)*two/sqn
               grady(i,j) = y(i,j)*(-one+x(i,j)**2+y(i,j)**2)
               grady(i,j) = grady(i,j)*two/sqn
               gradax(i,j) = zero
               graday(i,j) = zero
  130       continue
  140    continue

c        Kinetic Energy Part, Interior Points

         do 160 i = 2, nx
            do 150 j = 2, ny
               fkinx1 = (two/(hx*hx*sqn))*(x(i+1,j)-
     +                  x(i,j)*cos(hx*vpotx(i,j))+
     +                  y(i,j)*sin(hx*vpotx(i,j)))
               fkinx2 = (two/(hx*hx*sqn))*(y(i+1,j)-
     +                  y(i,j)*cos(hx*vpotx(i,j))-
     +                  x(i,j)*sin(hx*vpotx(i,j)))
               fkiny1 = (two/(hy*hy*sqn))*(x(i,j+1)-
     +                  x(i,j)*cos(hy*vpoty(i,j))+
     +                  y(i,j)*sin(hy*vpoty(i,j)))
               fkiny2 = (two/(hy*hy*sqn))*(y(i,j+1)-
     +                  y(i,j)*cos(hy*vpoty(i,j))-
     +                  x(i,j)*sin(hy*vpoty(i,j)))
               ffield = (vpotx(i,j)-vpotx(i,j+1))/hy +
     +                  (vpoty(i+1,j)-vpoty(i,j))/hx
               ffield = (two*(tkappa**2)/sqn)*ffield
               gradx(i,j) = gradx(i,j) - cos(hx*vpotx(i,j))*fkinx1 -
     +                      sin(hx*vpotx(i,j))*fkinx2 -
     +                      cos(hy*vpoty(i,j))*fkiny1 -
     +                      sin(hy*vpoty(i,j))*fkiny2
               grady(i,j) = grady(i,j) + sin(hx*vpotx(i,j))*fkinx1 -
     +                      cos(hx*vpotx(i,j))*fkinx2 +
     +                      sin(hy*vpoty(i,j))*fkiny1 -
     +                      cos(hy*vpoty(i,j))*fkiny2
               gradax(i,j) = gradax(i,j) + ffield/hy +
     +                       fkinx1*(hx*x(i,j)*sin(hx*vpotx(i,j))+
     +                       hx*y(i,j)*cos(hx*vpotx(i,j))) +
     +                       fkinx2*(hx*y(i,j)*sin(hx*vpotx(i,j))-
     +                       hx*x(i,j)*cos(hx*vpotx(i,j)))
               graday(i,j) = graday(i,j) - ffield/hx +
     +                       fkiny1*(hy*x(i,j)*sin(hy*vpoty(i,j))+
     +                       hy*y(i,j)*cos(hy*vpoty(i,j))) +
     +                       fkiny2*(hy*y(i,j)*sin(hy*vpoty(i,j))-
     +                       hy*x(i,j)*cos(hy*vpoty(i,j)))
               fkinx1 = (two/(hx*hx*sqn))*(x(i,j)-
     +                  x(i-1,j)*cos(hx*vpotx(i-1,j))+
     +                  y(i-1,j)*sin(hx*vpotx(i-1,j)))
               fkinx2 = (two/(hx*hx*sqn))*(y(i,j)-
     +                  y(i-1,j)*cos(hx*vpotx(i-1,j))-
     +                  x(i-1,j)*sin(hx*vpotx(i-1,j)))
               fkiny1 = (two/(hy*hy*sqn))*(x(i,j)-
     +                  x(i,j-1)*cos(hy*vpoty(i,j-1))+
     +                  y(i,j-1)*sin(hy*vpoty(i,j-1)))
               fkiny2 = (two/(hy*hy*sqn))*(y(i,j)-
     +                  y(i,j-1)*cos(hy*vpoty(i,j-1))-
     +                  x(i,j-1)*sin(hy*vpoty(i,j-1)))
               gradx(i,j) = gradx(i,j) + fkinx1 + fkiny1
               grady(i,j) = grady(i,j) + fkinx2 + fkiny2
               ffield = (vpotx(i,j-1)-vpotx(i,j))/hy +
     +                  (vpoty(i+1,j-1)-vpoty(i,j-1))/hx
               ffield = (two*(tkappa**2)/sqn)*ffield
               gradax(i,j) = gradax(i,j) - ffield/hy
               ffield = (vpotx(i-1,j)-vpotx(i-1,j+1))/hy +
     +                  (vpoty(i,j)-vpoty(i-1,j))/hx
               ffield = (two*(tkappa**2)/sqn)*ffield
               graday(i,j) = graday(i,j) + ffield/hx
  150       continue
  160    continue

c        Kinetic Energy Part, Boundary Points.

c        Bottom J = 1

         do 170 i = 2, nx
            fkinx1 = (two/(hx*hx*sqn))*(x(i+1,1)-
     +               x(i,1)*cos(hx*vpotx(i,1))+
     +               y(i,1)*sin(hx*vpotx(i,1)))
            fkinx2 = (two/(hx*hx*sqn))*(y(i+1,1)-
     +               y(i,1)*cos(hx*vpotx(i,1))-
     +               x(i,1)*sin(hx*vpotx(i,1)))
            fkiny1 = (two/(hy*hy*sqn))*(x(i,2)-
     +               x(i,1)*cos(hy*vpoty(i,1))+
     +               y(i,1)*sin(hy*vpoty(i,1)))
            fkiny2 = (two/(hy*hy*sqn))*(y(i,2)-
     +               y(i,1)*cos(hy*vpoty(i,1))-
     +               x(i,1)*sin(hy*vpoty(i,1)))
            ffield = (vpotx(i,1)-vpotx(i,2))/hy +
     +               (vpoty(i+1,1)-vpoty(i,1))/hx
            ffield = (two*(tkappa**2)/sqn)*ffield
            gradx(i,1) = gradx(i,1) - cos(hx*vpotx(i,1))*fkinx1 -
     +                   sin(hx*vpotx(i,1))*fkinx2 -
     +                   cos(hy*vpoty(i,1))*fkiny1 -
     +                   sin(hy*vpoty(i,1))*fkiny2
            grady(i,1) = grady(i,1) + sin(hx*vpotx(i,1))*fkinx1 -
     +                   cos(hx*vpotx(i,1))*fkinx2 +
     +                   sin(hy*vpoty(i,1))*fkiny1 -
     +                   cos(hy*vpoty(i,1))*fkiny2
            gradax(i,1) = gradax(i,1) + ffield/hy +
     +                    fkinx1*(hx*x(i,1)*sin(hx*vpotx(i,1))+
     +                    hx*y(i,1)*cos(hx*vpotx(i,1))) +
     +                    fkinx2*(hx*y(i,1)*sin(hx*vpotx(i,1))-
     +                    hx*x(i,1)*cos(hx*vpotx(i,1)))
            graday(i,1) = graday(i,1) - ffield/hx +
     +                    fkiny1*(hy*x(i,1)*sin(hy*vpoty(i,1))+
     +                    hy*y(i,1)*cos(hy*vpoty(i,1))) +
     +                    fkiny2*(hy*y(i,1)*sin(hy*vpoty(i,1))-
     +                    hy*x(i,1)*cos(hy*vpoty(i,1)))
            fkinx1 = (two/(hx*hx*sqn))*(x(i,1)-
     +               x(i-1,1)*cos(hx*vpotx(i-1,1))+
     +               y(i-1,1)*sin(hx*vpotx(i-1,1)))
            fkinx2 = (two/(hx*hx*sqn))*(y(i,1)-
     +               y(i-1,1)*cos(hx*vpotx(i-1,1))-
     +               x(i-1,1)*sin(hx*vpotx(i-1,1)))
            fkiny1 = (two/(hy*hy*sqn))*(x(i,ny+1)-
     +               x(i,ny)*cos(hy*vpoty(i,ny))+
     +               y(i,ny)*sin(hy*vpoty(i,ny)))
            fkiny2 = (two/(hy*hy*sqn))*(y(i,ny+1)-
     +               y(i,ny)*cos(hy*vpoty(i,ny))-
     +               x(i,ny)*sin(hy*vpoty(i,ny)))
            gradx(i,1) = gradx(i,1) + fkinx1 + fkiny1
            grady(i,1) = grady(i,1) + fkinx2 + fkiny2
            ffield = (vpotx(i,ny)-vpotx(i,ny+1))/hy +
     +               (vpoty(i+1,ny)-vpoty(i,ny))/hx
            ffield = (two*(tkappa**2)/sqn)*ffield
            gradax(i,1) = gradax(i,1) - ffield/hy
            ffield = (vpotx(i-1,1)-vpotx(i-1,2))/hy +
     +               (vpoty(i,1)-vpoty(i-1,1))/hx
            ffield = (two*(tkappa**2)/sqn)*ffield
            graday(i,1) = graday(i,1) + ffield/hx
  170    continue

c        Left I = 1.

         do 180 j = 2, ny
            fkinx1 = (two/(hx*hx*sqn))*(x(2,j)-
     +               x(1,j)*cos(hx*vpotx(1,j))+
     +               y(1,j)*sin(hx*vpotx(1,j)))
            fkinx2 = (two/(hx*hx*sqn))*(y(2,j)-
     +               y(1,j)*cos(hx*vpotx(1,j))-
     +               x(1,j)*sin(hx*vpotx(1,j)))
            fkiny1 = (two/(hy*hy*sqn))*(x(1,j+1)-
     +               x(1,j)*cos(hy*vpoty(1,j))+
     +               y(1,j)*sin(hy*vpoty(1,j)))
            fkiny2 = (two/(hy*hy*sqn))*(y(1,j+1)-
     +               y(1,j)*cos(hy*vpoty(1,j))-
     +               x(1,j)*sin(hy*vpoty(1,j)))
            ffield = (vpotx(1,j)-vpotx(1,j+1))/hy +
     +               (vpoty(2,j)-vpoty(1,j))/hx
            ffield = (two*(tkappa**2)/sqn)*ffield
            gradx(1,j) = gradx(1,j) - cos(hx*vpotx(1,j))*fkinx1 -
     +                   sin(hx*vpotx(1,j))*fkinx2 -
     +                   cos(hy*vpoty(1,j))*fkiny1 -
     +                   sin(hy*vpoty(1,j))*fkiny2
            grady(1,j) = grady(1,j) + sin(hx*vpotx(1,j))*fkinx1 -
     +                   cos(hx*vpotx(1,j))*fkinx2 +
     +                   sin(hy*vpoty(1,j))*fkiny1 -
     +                   cos(hy*vpoty(1,j))*fkiny2
            gradax(1,j) = gradax(1,j) + ffield/hy +
     +                    fkinx1*(hx*x(1,j)*sin(hx*vpotx(1,j))+
     +                    hx*y(1,j)*cos(hx*vpotx(1,j))) +
     +                    fkinx2*(hx*y(1,j)*sin(hx*vpotx(1,j))-
     +                    hx*x(1,j)*cos(hx*vpotx(1,j)))
            graday(1,j) = graday(1,j) - ffield/hx +
     +                    fkiny1*(hy*x(1,j)*sin(hy*vpoty(1,j))+
     +                    hy*y(1,j)*cos(hy*vpoty(1,j))) +
     +                    fkiny2*(hy*y(1,j)*sin(hy*vpoty(1,j))-
     +                    hy*x(1,j)*cos(hy*vpoty(1,j)))
            fkinx1 = (two/(hx*hx*sqn))*(x(nx+1,j)-
     +               x(nx,j)*cos(hx*vpotx(nx,j))+
     +               y(nx,j)*sin(hx*vpotx(nx,j)))
            fkinx2 = (two/(hx*hx*sqn))*(y(nx+1,j)-
     +               y(nx,j)*cos(hx*vpotx(nx,j))-
     +               x(nx,j)*sin(hx*vpotx(nx,j)))
            fkiny1 = (two/(hy*hy*sqn))*(x(1,j)-
     +               x(1,j-1)*cos(hy*vpoty(1,j-1))+
     +               y(1,j-1)*sin(hy*vpoty(1,j-1)))
            fkiny2 = (two/(hy*hy*sqn))*(y(1,j)-
     +               y(1,j-1)*cos(hy*vpoty(1,j-1))-
     +               x(1,j-1)*sin(hy*vpoty(1,j-1)))
            sfac = sin(two*pi*vornum*(j-one)/dble(ny))
            cfac = cos(two*pi*vornum*(j-one)/dble(ny))
            gradx(1,j) = gradx(1,j) + cfac*fkinx1 + sfac*fkinx2 + fkiny1
            grady(1,j) = grady(1,j) - sfac*fkinx1 + cfac*fkinx2 + fkiny2
            ffield = (vpotx(1,j-1)-vpotx(1,j))/hy +
     +               (vpoty(2,j-1)-vpoty(1,j-1))/hx
            ffield = (two*(tkappa**2)/sqn)*ffield
            gradax(1,j) = gradax(1,j) - ffield/hy
            ffield = (vpotx(nx,j)-vpotx(nx,j+1))/hy +
     +               (vpoty(nx+1,j)-vpoty(nx,j))/hx
            ffield = (two*(tkappa**2)/sqn)*ffield
            graday(1,j) = graday(1,j) + ffield/hx
  180    continue

c        Kinetic Energy Part, at origin (only needed in zero field).

         fkinx1 = (two/(hx*hx*sqn))*(x(2,1)-x(1,1)*cos(hx*vpotx(1,1))+
     +            y(1,1)*sin(hx*vpotx(1,1)))
         fkinx2 = (two/(hx*hx*sqn))*(y(2,1)-y(1,1)*cos(hx*vpotx(1,1))-
     +            x(1,1)*sin(hx*vpotx(1,1)))
         fkiny1 = (two/(hy*hy*sqn))*(x(1,2)-x(1,1)*cos(hy*vpoty(1,1))+
     +            y(1,1)*sin(hy*vpoty(1,1)))
         fkiny2 = (two/(hy*hy*sqn))*(y(1,2)-y(1,1)*cos(hy*vpoty(1,1))-
     +            x(1,1)*sin(hy*vpoty(1,1)))
         ffield = (vpotx(1,1)-vpotx(1,2))/hy +
     +            (vpoty(2,1)-vpoty(1,1))/hx
         ffield = (two*(tkappa**2)/sqn)*ffield
         gradx(1,1) = gradx(1,1) - cos(hx*vpotx(1,1))*fkinx1 -
     +                sin(hx*vpotx(1,1))*fkinx2 -
     +                cos(hy*vpoty(1,1))*fkiny1 -
     +                sin(hy*vpoty(1,1))*fkiny2
         grady(1,1) = grady(1,1) + sin(hx*vpotx(1,1))*fkinx1 -
     +                cos(hx*vpotx(1,1))*fkinx2 +
     +                sin(hy*vpoty(1,1))*fkiny1 -
     +                cos(hy*vpoty(1,1))*fkiny2
         gradax(1,1) = gradax(1,1) + ffield/hy +
     +                 fkinx1*(hx*x(1,1)*sin(hx*vpotx(1,1))+
     +                 hx*y(1,1)*cos(hx*vpotx(1,1))) +
     +                 fkinx2*(hx*y(1,1)*sin(hx*vpotx(1,1))-
     +                 hx*x(1,1)*cos(hx*vpotx(1,1)))
         graday(1,1) = graday(1,1) - ffield/hx +
     +                 fkiny1*(hy*x(1,1)*sin(hy*vpoty(1,1))+
     +                 hy*y(1,1)*cos(hy*vpoty(1,1))) +
     +                 fkiny2*(hy*y(1,1)*sin(hy*vpoty(1,1))-
     +                 hy*x(1,1)*cos(hy*vpoty(1,1)))
         fkinx1 = (two/(hx*hx*sqn))*(x(nx+1,1)-
     +            x(nx,1)*cos(hx*vpotx(nx,1))+
     +            y(nx,1)*sin(hx*vpotx(nx,1)))
         fkinx2 = (two/(hx*hx*sqn))*(y(nx+1,1)-
     +            y(nx,1)*cos(hx*vpotx(nx,1))-
     +            x(nx,1)*sin(hx*vpotx(nx,1)))
         fkiny1 = (two/(hy*hy*sqn))*(x(1,ny+1)-
     +            x(1,ny)*cos(hy*vpoty(1,ny))+
     +            y(1,ny)*sin(hy*vpoty(1,ny)))
         fkiny2 = (two/(hy*hy*sqn))*(y(1,ny+1)-
     +            y(1,ny)*cos(hy*vpoty(1,ny))-
     +            x(1,ny)*sin(hy*vpoty(1,ny)))
         gradx(1,1) = gradx(1,1) + fkinx1 + fkiny1
         grady(1,1) = grady(1,1) + fkinx2 + fkiny2
         ffield = (vpotx(1,ny)-vpotx(1,ny+1))/hy +
     +            (vpoty(2,ny)-vpoty(1,ny))/hx
         ffield = (two*(tkappa**2)/sqn)*ffield
         gradax(1,1) = gradax(1,1) - ffield/hy
         ffield = (vpotx(nx,1)-vpotx(nx,2))/hy +
     +            (vpoty(nx+1,1)-vpoty(nx,1))/hx
         ffield = (two*(tkappa**2)/sqn)*ffield
         graday(1,1) = graday(1,1) + ffield/hx
      end if

      end
