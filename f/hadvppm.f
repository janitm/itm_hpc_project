      subroutine hadvppm(nn,dt,dx,con,vel,mscl,flxarr,flux1,
     &                   flux2,saflux,fc1,fc2)


      include "camx.prm"

      real con(nn),vel(nn),flxarr(nn),mscl(nn),saflux(nn)
      real*8 flux1,flux2

      real fc1(mx1d),fc2(mx1d)
      parameter (TWO3RDS=2./3.)

      do i = 1,nn
        fm(i) = 0.
        fp(i) = 0.
        fc1(i) = 0.
        fc2(i) = 0.
      enddo

      cm(2)  = con(2)
      cm(nn) = con(nn-1)
      cm(3)    = (con(3) + con(2))/2.
      cm(nn-1) = (con(nn-1) + con(nn-2))/2.


      do i = 3,nn-2
        dc(i) = 0.5*(con(i+1) - con(i-1))
      
        if ((con(i+1) - con(i))*(con(i) - con(i-1)).gt.0.) then
          dc(i) = sign(1.,dc(i))*min(
     &                               abs(dc(i)),
     &                               2.*abs(con(i+1) - con(i)),
     &                               2.*abs(con(i) - con(i-1)))
        else
          dc(i) = 0.
        endif
      enddo

      do i = 3,nn-3
        cm(i+1) = con(i) + 
     &            0.5*(con(i+1) - con(i)) + (dc(i) - dc(i+1))/6.
      enddo

      do i = 2,nn-1
        cr(i) = cm(i+1)
        cl(i) = cm(i)
      enddo

      do i = 2,nn-1
        if ((cr(i) - con(i))*(con(i) - cl(i)).gt.0.) then
          dc(i) = cr(i) - cl(i)
          c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))

          if (dc(i)*c6(i) .gt. dc(i)*dc(i)) then
            cl(i) = 3.*con(i) - 2.*cr(i)
          elseif (-dc(i)*dc(i) .gt. dc(i)*c6(i)) then
            cr(i) = 3.*con(i) - 2.*cl(i)
          endif
        else
          cl(i) = con(i)
          cr(i) = con(i)
        endif
        dc(i) = cr(i) - cl(i)
        c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
      enddo

      do i = 2,nn-1
        x = max(0., -vel(i-1)*(dt/dx))
        fm(i) = x*(cl(i) + 0.5*x*(dc(i) + c6(i)*(1. - TWO3RDS*x)))
        if (x.ge.1) write(*,*)'Courant number is bigger than 1',x ! jgj 10/6/06
        x = max(0., vel(i)*(dt/dx))
        if (x.ge.1) write(*,*)'Courant number is bigger than 1',x ! jgj 10/6/06
        fp(i) = x*(cr(i) - 0.5*x*(dc(i) - c6(i)*(1. - TWO3RDS*x)))
      enddo

      if (vel(1).gt.0.) then
        x = vel(1)*(dt/dx)
        fp(1) = x*con(1)
      endif

      if (vel(nn-1).lt.0.) then
        x = -vel(nn-1)*(dt/dx)
        fm(nn) = x*con(nn)
      endif

      flxarr(1) = (fp(1) - fm(2))*(dx/dt)
      saflux(1) = flxarr(1)*(dt/dx)
      do i = 2,nn-1
        flxarr(i) = (fp(i) - fm(i+1))*(dx/dt)
        con(i) = con(i) - mscl(i)*(flxarr(i) - flxarr(i-1))*(dt/dx)
        saflux(i) = flxarr(i)*(dt/dx)

        fc1(i) =   mscl(i)*flxarr(i-1)*(dt/dx)
        fc2(i) = - mscl(i)*flxarr(i)*(dt/dx)
      enddo

      flux1 = mscl(2)*flxarr(1)
      flux2 = mscl(nn-1)*flxarr(nn-1)

      return
      end
