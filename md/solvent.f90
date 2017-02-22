

module solvent

	implicit none

        integer npts,ncpu
        double precision w(nmax),gr2(2,nmax)
        double precision alp(nmax),g0(nmax),x0(nmax)
        double precision As(nmax),Bs(nmax),vtot(nmax)

	contains
	
	! set parameters for the solvent model
        subroutine setsolvent()
        integer nmax
        parameter(nmax=10000)
C solvent values
        integer npts,ncpu
        double precision w(nmax),gr2(2,nmax)
        double precision alp(nmax),g0(nmax),x0(nmax)
        double precision As(nmax),Bs(nmax),vtot(nmax)
        double precision iin(2)
C
        double precision gr(3,nmax)
        integer i,j
C
        common /solvent/ gr2,alp,x0,g0,As,Bs,vtot,w,npts,ncpu
C
        npts=1000
C
        As(1)=7.85890042d5*12.d0
        As(2)=6.91773368d4*12.d0
C
        Bs(1)=6.36687196d2*6.d0
        Bs(2)=1.16264660d2*6.d0
C
        x0(1)= 4.0665871952666128d0     
        g0(1)= 1.1776140774892252d0     
        alp(1)=0.6222126557254711d0     
        x0(2)= 3.0531410266640604d0     
        g0(2)= 1.3845336480728017d0     
        alp(2)=2.3812244519571273d0     
        iin(1)=1
        iin(2)=1
C
        gr(1,1)=x0(1)-dsqrt(g0(1)/alp(1))
        gr(2,1)=x0(1)+dble(iin(1))*dsqrt((g0(1)-1.d0)/alp(1))
        gr(3,1)=x0(1)+dsqrt(g0(1)/alp(1))
        gr(1,2)=x0(2)-dsqrt(g0(2)/alp(2))
        gr(2,2)=x0(2)+dble(iin(2))*dsqrt((g0(2)-1.d0)/alp(2))
        gr(3,2)=x0(2)+dsqrt(g0(2)/alp(2))
C
        do i=1,2
          do j=1,2
            gr2(j,i)=gr(j,i)**2
          enddo
        enddo
C
        do i=1,2
          w(i)=(gr(3,i)-gr(1,i))/2.d0
        enddo
C
        do i=1,2
          vtot(i)=16.d0/3.d0*3.1415927d0*w(i)*g0(i)/dble(npts)*0.0334d0
        enddo
C
      return
      endsubroutine setsolvent

endmodule solvent

