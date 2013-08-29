      program runme
c
c    Wrapper for running hadvppm for use in HPC project
c    Ben Murphy (August, 2013)
c
c
      integer nn
      parameter(nn=90)  !Length of Dimension
      real c1d(nn), v1d(nn), cinit(nn), m1d(nn), dt, dx
      real flxarr(nn), flux1, flux2
      integer ii


      !Set up input scalars, vectors, and arrays
      do ii = 1,nn
        c1d(ii) = 1.0 !Concentration (micrograms / m3)
	cinit(ii) = c1d(ii)
	v1d(ii) = 0.5 !Wind Speed (m/s)
	m1d(ii) = 1.0 !Map Scale Factor (important for projecting
	              !coordinates onto grid, set 1.0 for this 
		      !exercise
        flxarr(ii) = 0.0  !Output Vector. Set to 0 intitially
      enddo

      dt = 900     !Time Step (seconds)
      dx = 36000.0 !Grid size (m)
      flux1 = 0.0  !Output - flux at start boundary
      flux2 = 0.0  !Output - flux at end boundary

      !CALL 1-D ADVECTION SUBROUTINE!
      print *,'Calling hadvppm...'
      call hadvppm(nn,dt,dx,c1d,v1d,m1d,flxarr,flux1,
     &                    flux2)
      flxarr(nn) = 0.0  !The flux is BETWEEN grid cells, so for the last
                        !cell it equals 0.

      !Output all fo the data to a text file
      print *,'Writing output data...'
      open(unit=20,file='Flux_Out.dat')
      
      write(20,'(A4,2x,A9,2x,A9,2x,A4)'), 'Cell','Start Con','End Con','Flux'
      do ii = 1,nn
        write(20,'(I2,4x,E9.4,2x,E9.4,2x,E9.4)'), ii, cinit(ii), c1d(ii), flxarr(ii)
      enddo
      
      close(20)


      end
