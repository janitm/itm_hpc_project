
cdbg        subroutine mnratio(ii,ngas,naero,c1d,c2d,i1,i2,MX1D,mn,j,k,flag)
        subroutine mnratio(ii,c1d,c2d,i1,i2,mn,order,nspc,j,k,flag)

        include "camx.prm"
        include "chmstry.com"

cdbg        integer ii,ngas,naero
        integer ii
        integer i,spec
        integer j,k
        integer flag ! 1 => x-advection, 2 => y-advection
        integer order(nspc)
        real c1d(MX1D),c2d(MX1D),mn(MX1D)
        real eps

        parameter (eps = 1.0e-20)
c
        if (ii.gt.ngas) then
           spec=mod(ii-ngas,naero)
           if (spec.eq.1) then
             l = 0
             do i = i1,i2
               l = l + 1
               c2d(l)=c1d(l)
               if (c2d(l) .le. 0.0) then
                 if (c2d(l) .gt. -eps) then
                   c2d(l) = eps
                 else
                   print*,'In mnratio, number concentrations are less
     &                   than -eps'
                   print*,'Species=',spname(order(ii))
                   print*,'Concentration=',c2d(l)
                   print*,'Coordinate='
                   if (flag.eq.1) then
                     print*,(i,j,k)
                   elseif (flag.eq.2) then
                     print*,(j,i,k)
                   endif
                   pause
                 endif
               endif
               mn(l)=0.0
             enddo
           else
             l = 0
             do i = i1,i2
               l = l + 1
               mn(l)=c1d(l)/c2d(l)
             enddo
           endif
        endif

        return
        end

c 0000000000000000000000000000000000000000000000000000000000000000000


        subroutine hadvppm2(ii,ngas,naero,MX1D,mn,nflxarr,
     &                  nn,dt,dx,con,mscl,flxarr,flux1,flux2,
     &                  saflux,fc1,fc2)

        integer ii,ngas,naero
        integer jj,spec
        real con(nn),mn(nn),nflxarr(nn),flxarr(nn),mscl(nn),saflux(nn)
        real*8 flux1,flux2

         if (nflxarr(1).gt.0) then
            flxarr(1) = nflxarr(1)*mn(1)
         else
            flxarr(1) = nflxarr(1)*mn(2)
         endif

         do jj = 2,nn-1
            if (nflxarr(jj).gt.0) then
               flxarr(jj) = nflxarr(jj)*mn(jj)
            else
               flxarr(jj) = nflxarr(jj)*mn(jj+1)
            endif
            con(jj) = con(jj) - mscl(jj)*(flxarr(jj) - flxarr(jj-1))
     &              *(dt/dx)
            saflux(jj) = flxarr(jj)*(dt/dx)
         enddo

          do jj = 2,nn-1
             fc1(jj) =   mscl(jj)*flxarr(jj-1)*(dt/dx)
             fc2(jj) = - mscl(jj)*flxarr(jj)*(dt/dx)
          enddo

             flux1 = mscl(2)*flxarr(1)
             flux2 = mscl(nn-1)*flxarr(nn-1)

            return
            end
