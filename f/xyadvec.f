      subroutine xyadvec(igrd,xyordr,ncol,nrow,nlay,nspc,nsen,nadv,
     &                   deltat,dx,dy,windu,windv,depth,mapscl,conc,
     &                   fluxes,sens,tarray2,isaptr,species)


      include "camx.prm"
      include "bndary.com"
      include "chmstry.com"
      include "filunit.com"
      include "flags.com"


      integer xyordr
      integer nadv(nlay)
      dimension conc(ncol,nrow,nlay,nspc)
      dimension sens(ncol,nrow,nlay,nsen)
      real windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &     depth(ncol,nrow,nlay),
     &     mapscl(ncol,nrow),
     &     dx(nrow)
      real c1d(MX1D),v1d(MX1D),flxarr(MX1D),m1d(MX1D),saflux(MX1D)
      real fluxo3(MX1D), fluxvoc(MX1D), fluxnox(MX1D)
      real cnco3(MX1D), cncvoc(MX1D), cncnox(MX1D)
      real c1d0(MX1D),fpc(MX1D),fmc(MX1D)
      real sen1d(MX1D,MXTRSP)
      real*8 fluxes(nspc,11),flux1,flux2
      real*8 fluxtmp(MXSPEC,8,MXLAYA)
      real c2d(MX1D),mn(MX1D) ! by jgj 8/13/06
      dimension tarray2(2)
      character*10 species(nspc) ! by jgj 7/18/06

      integer ii,jj,kk ! for generating the order of an advection by jgj 8/8/06
      integer order(nspc) ! the order of species by jgj 8/8/06
      integer nsect ! number of size section by jgj 8/8/06
      integer spec ! variable for calculating species by jgj 8/14/06
      integer lmod
      real nflxarr(MX1D) ! flxarr for number species by jgj 8/14/06

      nsect = 43
      do kk=1,ngas
         order(kk) = kk
      enddo

      kk=ngas

      do ii=1,nsect
         do jj=1,naero
            kk=kk+1
            if (jj.eq.1) then
               order(kk) = ngas + (naero-1)*nsect + ii
            else
               order(kk) = ngas + (jj-2)*nsect + ii
            endif
         enddo
      enddo


      do i=1,nspc
        do j=1,8
          do k=1,nlay
            fluxtmp(i,j,k) = 0.
          enddo
        enddo
      enddo

      if( xyordr .eq. 0 ) goto 200

 100  continue

      write(*,'(a20,$)') 'x advection ......'
      write(iout,'(a20,$)') 'x advection ......'
      tcpu = dtime(tarray2)

      do 20 k = 1,nlay
        dtuse = deltat/nadv(k)
        do 21 istep = 1,nadv(k)

        do 10 j = 2,nrow-1
          i1 = 1
          i2 = ncol
          if (igrd.eq.1) then
            if (ibeg(j).eq.-999) goto 10
            i1 = ibeg(j) - 1
            i2 = iend(j) + 1
          endif
          l = 0
          do i = i1,i2-1
            l = l + 1
            v1d(l) = 2.*windu(i,j,k)/(mapscl(i+1,j) + mapscl(i,j))
            m1d(l) = mapscl(i,j)*mapscl(i,j)
          enddo
          l = l + 1
          v1d(l) = windu(i2,j,k)
          m1d(l) = mapscl(i2,j)*mapscl(i2,j)
          nn = i2 - i1 + 1

          do 22 ii = 1,nspc
            ispc = order(ii)

            l = 0
            do i = i1,i2
              l = l + 1
              c1d(l) = conc(i,j,k,ispc)*dy*depth(i,j,k)
            enddo
            if (igrd.eq.1 .and. v1d(1).lt.0.) c1d(1) = c1d(2)
            if (igrd.eq.1 .and. v1d(nn-1).gt.0.) c1d(nn) = c1d(nn-1)

            call mnratio(ii,c1d,c2d,i1,i2,mn,order,nspc,j,k,1) ! by jgj 8/13/06

            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then ! An number concentration
                  do jj = 1,nn-1
                     nflxarr(jj) = 0
                  enddo
               else
                  call hadvppm2(ii,ngas,naero,MX1D,mn,nflxarr,
     &             nn,dtuse,dx(j),c1d,m1d,flxarr,flux1,flux2,
     &             saflux,fc1,fc2)
                  goto 400
               endif
            endif

          CALL 1-D ADVECTION SUBROUTINE!
            call hadvppm(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                           flux2,saflux,fc1,fc2)

            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then
                  do jj = 1,nn-1
                     nflxarr(jj) = flxarr(jj)
                  enddo
               endif
            endif

 400        continue

            l = 1
            do i = i1+1,i2-1
              l = l + 1
              conc(i,j,k,ispc) = c1d(l)/dy/depth(i,j,k)
            enddo
  22      continue

  10    continue

  21  continue
  20  continue

      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)

      if (xyordr.eq.0) goto 300

 200  continue

      write(*,'(a20,$)') 'y advection ......'
      write(iout,'(a20,$)') 'y advection ......'

      do 40 k = 1,nlay
        dtuse = deltat/nadv(k)
        do 41 istep = 1,nadv(k)

        do 30 i = 2,ncol-1
          j1 = 1
          j2 = nrow
          if (igrd.eq.1) then
            if (jbeg(i).eq.-999) goto 30
            j1 = jbeg(i) - 1
            j2 = jend(i) + 1
          endif
          l = 0
          do j = j1,j2-1
            l = l + 1
            v1d(l) = 2.*windv(i,j,k)/(mapscl(i,j+1) + mapscl(i,j))
            m1d(l) = mapscl(i,j)*mapscl(i,j)
          enddo
          l = l + 1
          v1d(l) = windv(i,j2,k)
          m1d(l) = mapscl(i,j2)*mapscl(i,j2)
          nn = j2 - j1 + 1

          do 42 ii = 1,nspc
            ispc = order(ii)

            l = 0
            do j = j1,j2
              l = l + 1

              c1d(l) = conc(i,j,k,ispc)*dx(j)*depth(i,j,k)
            enddo
            if (igrd.eq.1 .and. v1d(1).lt.0.) c1d(1) = c1d(2)
            if (igrd.eq.1 .and. v1d(nn-1).gt.0.) c1d(nn) = c1d(nn-1)

            call mnratio(ii,c1d,c2d,j1,j2,mn,order,nspc,i,k,2) ! by jgj 8/13/06

            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then
                  do jj = 1,nn-1
                     nflxarr(jj) = 0
                  enddo
               else
                  call hadvppm2(ii,ngas,naero,MX1D,mn,nflxarr,
     &             nn,dtuse,dy,c1d,m1d,flxarr,flux1,flux2,
     &             saflux,fc1,fc2)
                  goto 500
               endif
            endif

            if (iadvct.eq.2) then
              call hadvbot(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                 flux2,saflux,fpc,fmc,fc1,fc2)
            elseif( iadvct .eq. 3) then
              call hadvppm(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                         flux2,saflux,fc1,fc2)
            endif

            if (ii.gt.ngas) then
               spec=mod(ii-ngas,naero)
               if (spec.eq.1) then
                  do jj = 1,nn-1
                     nflxarr(jj) = flxarr(jj)
                  enddo
               endif
            endif

 500        continue

            l = 1
            do j = j1+1,j2-1
              l = l+1
              conc(i,j,k,ispc) = c1d(l)/dx(j)/depth(i,j,k)
            enddo

  42      continue

  30    continue

  41  continue
  40  continue

      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)

      if (xyordr.eq.0) goto 100

 300  continue

      call flush(6)
      call flush(iout)
      return
      end
