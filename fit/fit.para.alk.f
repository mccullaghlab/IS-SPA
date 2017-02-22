CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM fit_gr
        integer nmax,mpts
        parameter(nmax=100,mpts=5000000)
        double precision y(mpts),x(3,mpts)
        double precision xpos(3,nmax)
        double precision r12(2,nmax),r22(2,nmax)
        double precision r12p(2,nmax),r22p(2,nmax)
        double precision r10(2,nmax),r20(2,nmax)
        double precision r10p(2,nmax),r20p(2,nmax)
        double precision a(nmax,nmax),b(nmax),dy(nmax)
        double precision val(nmax),sval(nmax),valp(nmax)
        double precision chi2,chi2p,lam,fnow,r0
        double precision gnow,x1,x2,x3,r2
        integer ityp(nmax)
        integer i,j,k,igo,nval,npts,npar,ntyp,jt,jgo,idx
        double precision dx,dy1,dy2,dy3,sv1,sv2,sv0,v1,v2,v0
        double precision alp(nmax),alpp(nmax),dalp(3,nmax),salp(nmax)
C define initial values to the fit parameters
        ntyp=2
        nval=ntyp*3
        val(1)=2.8789d0
        val(2)=4.04137d0
        val(3)=3.65d0
        val(4)=2.37789d0
        val(5)=3.50257d0
        val(6)=3.15d0
C read in the molecule coordinates and types
        open(20,FILE='XXX',STATUS='old')
        read(20,*) npar
        do i=1,npar
          read(20,*) ityp(i),xpos(1,i),xpos(2,i),xpos(3,i)
        enddo
        close(20)
C read in the distribution function
        open(20,FILE='YYY',STATUS='old')
        read(20,*) npts
        do i=1,npts
          read(20,*) x(1,i),x(2,i),x(3,i),y(i)
        enddo
C loop over epsilon values
        dx=0.05d0*dsqrt(3.d0)
        do idx=1,4
          dx=dx*0.5d0
C define the distance parameters given the fitted parameters
        do i=1,ntyp
          r10(1,i)=val(3*i-2)-dx
          r10(2,i)=val(3*i-2)+dx
          r10(1,ntyp+i)=val(3*i-2)-dx
          r10(2,ntyp+i)=val(3*i-2)+dx
        enddo
        do i=1,ntyp
          r20(1,i)=val(3*i-1)-dx
          r20(2,i)=val(3*i-1)+dx
          r20(1,ntyp+i)=2.d0*val(3*i)-val(3*i-2)
          r20(1,ntyp+i)=r20(1,ntyp+i)-dx
          r20(2,ntyp+i)=r20(1,ntyp+i)+2.d0*dx
        enddo
        do i=1,2*ntyp
          r12(1,i)=r10(1,i)**2
          r12(2,i)=r10(2,i)**2
          r22(1,i)=r20(1,i)**2
          r22(2,i)=r20(2,i)**2
        enddo
        do i=1,ntyp
          v1=val(i*3-2)
          v2=val(i*3-1)
          v0=val(i*3)
          alp(i)=1.d0/((v1-v0)*(v1-v0)-(v2-v0)*(v2-v0))
          dalp(1,i)=-2.d0*(v1-v0)*alp(i)*alp(i)
          dalp(2,i)= 2.d0*(v2-v0)*alp(i)*alp(i)
          dalp(3,i)= 2.d0*(v1-v2)*alp(i)*alp(i)
        enddo
C calculate initial chi2
        chi2=0.d0
        do i=1,npts
          j=1
          jgo=1
          gnow=1.d0
          do while (jgo.eq.1)
            jt=ityp(j)
            if(j.eq.npar) jgo=0
            x1=xpos(1,j)-x(1,i)
            x2=xpos(2,j)-x(2,i)
            x3=xpos(3,j)-x(3,i)
            r2=x1*x1+x2*x2+x3*x3
            if(r2.lt.r12(1,jt)) then
              jgo=-1
            elseif(r2.lt.r22(2,jt)) then
              v1=val(jt*3-2)
              v2=val(jt*3-1)
              v0=val(jt*3)
              r0=dsqrt(r2)
              if(r2.gt.r22(1,jt)) then
                gnow=gnow*(1.d0+alp(jt)/4.d0/dx*(r0-v2-dx)**2*
     &               (r0-2.d0*v0+v2))
              elseif(r2.lt.r12(2,jt)) then
                gnow=-gnow*alp(jt)/4.d0/dx*(r0-v1+dx)**2*
     &               (r0-2.d0*v0+v1)
              else
                gnow=gnow*alp(jt)*((v1-v0)**2-(r0-v0)**2)
              endif
            endif
            j=j+1
          enddo
          if(jgo.eq.0) then
            gnow=y(i)-gnow
            chi2=chi2+gnow*gnow
          else
            chi2=chi2+y(i)*y(i)
          endif
        enddo
C start the fit
        lam=0.01d0
        igo=2
        do while(igo.gt.0)
          do i=1,nval
            b(i)=0.d0
            do j=1,nval
              a(j,i)=0.d0
            enddo
          enddo
C populate the a and b matrices
          do i=1,npts
            do j=1,nval
              dy(j)=0.d0
            enddo
            gnow=1.d0
            j=1
            jgo=1
            do while (jgo.eq.1)
              jt=ityp(j)
              if(j.eq.npar) jgo=0
              x1=xpos(1,j)-x(1,i)
              x2=xpos(2,j)-x(2,i)
              x3=xpos(3,j)-x(3,i)
              r2=x1*x1+x2*x2+x3*x3
              if(r2.lt.r12(1,jt)) then
                jgo=-1
              elseif(r2.lt.r22(2,jt)) then
                v1=val(jt*3-2)
                v2=val(jt*3-1)
                v0=val(jt*3)
                r0=dsqrt(r2)
                if(r2.gt.r22(1,jt)) then
                  fnow=alp(jt)/4.d0/dx*(r0-v2-dx)**2*
     &                 (r0-2.d0*v0+v2)
                  dy1=dalp(1,jt)/alp(jt)
                  dy2=dalp(2,jt)/alp(jt)
     &                 -2.d0/(r0-v2-dx)+1.d0/(r0-2.d0*v0+v2)
                  dy3=dalp(3,jt)/alp(jt)
     &                 -2.d0/(r0-2.d0*v0+v2)
                  dy1=dy1*fnow/(fnow+1.d0)
                  dy2=dy2*fnow/(fnow+1.d0)
                  dy3=dy3*fnow/(fnow+1.d0)
                  fnow=fnow+1.d0
                elseif(r2.lt.r12(2,jt)) then
                  fnow=-alp(jt)*(r0-v1+dx)**2*(r0-2.d0*v0+v1)/4.d0/dx
                  dy1=dalp(1,jt)/alp(jt)-2.d0/(r0-v1+dx)+
     &                 1.d0/(r0-2.d0*v0+v1)
                  dy2=dalp(2,jt)/alp(jt)
                  dy3=dalp(3,jt)/alp(jt)-2.d0/(r0-2.d0*v0+v1)
                else
                  fnow=alp(jt)*((v1-v0)**2-(r0-v0)**2)
                  dy1=dalp(1,jt)/alp(jt)+2.d0*(v1-v0)/
     &                 ((v1-v0)**2-(r0-v0)**2)
                  dy2=dalp(2,jt)/alp(jt)
                  dy3=dalp(3,jt)/alp(jt)+2.d0*(r0-v1)/
     &                 ((v1-v0)**2-(r0-v0)**2)
                endif
                gnow=gnow*fnow
                dy(jt*3-2)=dy(jt*3-2)+dy1
                dy(jt*3-1)=dy(jt*3-1)+dy2
                dy(jt*3)=dy(jt*3)+dy3
              endif
              j=j+1
            enddo
            if(jgo.eq.0) then
              do j=1,nval
                dy(j)=dy(j)*gnow
              enddo
              gnow=y(i)-gnow
              do j=1,nval
                b(j)=b(j)+dy(j)*gnow
                a(j,j)=a(j,j)+dy(j)*dy(j)*(1.d0+lam)
                do k=j+1,nval
                  a(j,k)=a(j,k)+dy(j)*dy(k)
                enddo
              enddo
            endif
          enddo
          do j=1,nval
            do k=j+1,nval
              a(k,j)=a(j,k)
            enddo
          enddo
C invert
          call inv(a,b,nval)
C populate temporary values
          do i=1,nval
            valp(i)=val(i)+b(i)
          enddo
          do i=1,ntyp
            r10p(1,i)=valp(3*i-2)-dx
            r10p(2,i)=valp(3*i-2)+dx
            r10p(1,ntyp+i)=valp(3*i-2)-dx
            r10p(2,ntyp+i)=valp(3*i-2)+dx
          enddo
          do i=1,ntyp
            r20p(1,i)=valp(3*i-1)-dx
            r20p(2,i)=valp(3*i-1)+dx
            r20p(1,ntyp+i)=2.d0*valp(3*i)-valp(3*i-2)
            r20p(1,ntyp+i)=r20p(1,ntyp+i)-dx
            r20p(2,ntyp+i)=r20p(1,ntyp+i)+2.d0*dx
          enddo
          do i=1,2*ntyp
            r12p(1,i)=r10p(1,i)**2
            r12p(2,i)=r10p(2,i)**2
            r22p(1,i)=r20p(1,i)**2
            r22p(2,i)=r20p(2,i)**2
          enddo
          do i=1,ntyp
            v1=valp(i*3-2)
            v2=valp(i*3-1)
            v0=valp(i*3)
          enddo
C calculate new chi2
          chi2p=0.d0
          do i=1,npts
            gnow=1.d0
            jgo=1
            j=1
            do while (jgo.eq.1)
              jt=ityp(j)
              if(j.eq.npar) jgo=0
              x1=xpos(1,j)-x(1,i)
              x2=xpos(2,j)-x(2,i)
              x3=xpos(3,j)-x(3,i)
              r2=x1*x1+x2*x2+x3*x3
              if(r2.lt.r12p(1,jt)) then
                jgo=-1
              elseif(r2.lt.r22p(2,jt)) then
                v1=valp(jt*3-2)
                v2=valp(jt*3-1)
                v0=valp(jt*3)
                r0=dsqrt(r2)
                if(r2.gt.r22p(1,jt)) then
                  gnow=gnow*(1.d0+alpp(jt)/4.d0/dx*(r0-v2-dx)**2*
     &                 (r0-2.d0*v0+v2))
                elseif(r2.lt.r12p(2,jt)) then
                  gnow=-gnow*alpp(jt)/4.d0/dx*(r0-v1+dx)**2*
     &                 (r0-2.d0*v0+v1)
                else
                  gnow=gnow*alpp(jt)*((v1-v0)**2-(r0-v0)**2)
                endif
              endif
              j=j+1
            enddo
            if(jgo.eq.0) then
              gnow=y(i)-gnow
              chi2p=chi2p+gnow*gnow
            else
              chi2p=chi2p+y(i)*y(i)
            endif
          enddo
C increment
          write(6,*) ''
          write(6,*) lam
          write(6,*) chi2,chi2p
          do i=1,nval
            write(6,*) val(i),valp(i)
          enddo
          write(6,*) ''
          if(chi2p.gt.chi2) then
            lam=lam*10.d0
            igo=2
          else
            lam=lam*0.1d0
            do i=1,nval
              val(i)=valp(i)
            enddo
            do i=1,2*ntyp
              r12(1,i)=r12p(1,i)
              r12(2,i)=r12p(2,i)
              r22(1,i)=r22p(1,i)
              r22(2,i)=r22p(2,i)
              r10(1,i)=r10p(1,i)
              r10(2,i)=r10p(2,i)
              r20(1,i)=r20p(1,i)
              r20(2,i)=r20p(2,i)
            enddo
            do i=1,ntyp
              alp(i)=alpp(i)
              v1=val(i*3-2)
              v2=val(i*3-1)
              v0=val(i*3)
              dalp(1,i)=-2.d0*(v1-v0)*alp(i)*alp(i)
              dalp(2,i)= 2.d0*(v2-v0)*alp(i)*alp(i)
              dalp(3,i)= 2.d0*(v1-v2)*alp(i)*alp(i)
            enddo
            if((chi2-chi2p)/chi2p.lt.1.d-5) igo=igo-1
            chi2=chi2p
          endif
        enddo
C
        write(6,*) 'chi2    =',chi2
        write(6,*) 'Cx1     =',val(1)
        write(6,*) 'Cx2     =',val(2)
        write(6,*) 'Cx0     =',val(3)
        write(6,*) 'Hx1     =',val(4)
        write(6,*) 'Hx2     =',val(5)
        write(6,*) 'Hx0     =',val(6)
        enddo
C
        lam=0.d0
        do i=1,nval
          b(i)=0.d0
          do j=1,nval
            a(j,i)=0.d0
          enddo
        enddo
C
        lam=0.d0
        dx=0.d0
        chi2=0.d0
        do i=1,ntyp
          r10(1,i)=val(3*i-2)-dx
          r10(2,i)=val(3*i-2)+dx
          r10(1,ntyp+i)=val(3*i-2)-dx
          r10(2,ntyp+i)=val(3*i-2)+dx
        enddo
        do i=1,ntyp
          r20(1,i)=val(3*i-1)-dx
          r20(2,i)=val(3*i-1)+dx
          r20(1,ntyp+i)=2.d0*val(3*i)-val(3*i-2)
          r20(1,ntyp+i)=r20(1,ntyp+i)-dx
          r20(2,ntyp+i)=r20(1,ntyp+i)+2.d0*dx
        enddo
        do i=1,2*ntyp
          r12(1,i)=r10(1,i)**2
          r12(2,i)=r10(2,i)**2
          r22(1,i)=r20(1,i)**2
          r22(2,i)=r20(2,i)**2
        enddo
        do i=1,npts
          do j=1,nval
            dy(j)=0.d0
          enddo
          gnow=1.d0
          j=1
          jgo=1
          do while (jgo.eq.1)
            jt=ityp(j)
            if(j.eq.npar) jgo=0
            x1=xpos(1,j)-x(1,i)
            x2=xpos(2,j)-x(2,i)
            x3=xpos(3,j)-x(3,i)
            r2=x1*x1+x2*x2+x3*x3
            if(r2.lt.r12(1,jt)) then
              jgo=-1
            elseif(r2.lt.r22(2,jt)) then
              v1=val(jt*3-2)
              v2=val(jt*3-1)
              v0=val(jt*3)
              r0=dsqrt(r2)
              if(r2.gt.r22(1,jt)) then
                fnow=alp(jt)/4.d0/dx*(r0-v2-dx)**2*
     &               (r0-2.d0*v0+v2)
                dy1=dalp(1,jt)/alp(jt)
                dy2=dalp(2,jt)/alp(jt)
     &               -2.d0/(r0-v2-dx)+1.d0/(r0-2.d0*v0+v2)
                dy3=dalp(3,jt)/alp(jt)
     &               -2.d0/(r0-2.d0*v0+v2)
                dy1=dy1*fnow/(fnow+1.d0)
                dy2=dy2*fnow/(fnow+1.d0)
                dy3=dy3*fnow/(fnow+1.d0)
                fnow=fnow+1.d0
              elseif(r2.lt.r12(2,jt)) then
                fnow=-alp(jt)*(r0-v1+dx)**2*(r0-2.d0*v0+v1)/4.d0/dx
                dy1=dalp(1,jt)/alp(jt)-2.d0/(r0-v1+dx)+
     &               1.d0/(r0-2.d0*v0+v1)
                dy2=dalp(2,jt)/alp(jt)
                dy3=dalp(3,jt)/alp(jt)-2.d0/(r0-2.d0*v0+v1)
              else
                fnow=alp(jt)*((v1-v0)**2-(r0-v0)**2)
                dy1=dalp(1,jt)/alp(jt)+2.d0*(v1-v0)/
     &               ((v1-v0)**2-(r0-v0)**2)
                dy2=dalp(2,jt)/alp(jt)
                dy3=dalp(3,jt)/alp(jt)+2.d0*(r0-v1)/
     &               ((v1-v0)**2-(r0-v0)**2)
              endif
              gnow=gnow*fnow
              dy(jt*3-2)=dy(jt*3-2)+dy1
              dy(jt*3-1)=dy(jt*3-1)+dy2
              dy(jt*3)=dy(jt*3)+dy3
            endif
            j=j+1
          enddo
          if(jgo.eq.0) then
            do j=1,nval
              dy(j)=dy(j)*gnow
            enddo
            gnow=y(i)-gnow
            do j=1,nval
              b(j)=b(j)+dy(j)*gnow
              a(j,j)=a(j,j)+dy(j)*dy(j)*(1.d0+lam)
              do k=j+1,nval
                a(j,k)=a(j,k)+dy(j)*dy(k)
              enddo
            enddo
            chi2=chi2+gnow*gnow
          else
            chi2=chi2+y(i)*y(i)
          endif
        enddo
        do j=1,nval
          do k=j+1,nval
            a(k,j)=a(j,k)
          enddo
        enddo
C invert
        call inv(a,b,nval)
        do i=1,nval
          sval(i)=dsqrt(a(i,i)*chi2/dble(npts-nval))
        enddo
        do i=1,ntyp
          v1=val(i*3-2)
          v2=val(i*3-1)
          v0=val(i*3)
          sv1=sval(i*3-2)**2
          sv2=sval(i*3-1)**2
          sv0=sval(i*3)**2
          salp(i)=2.d0*alp(i)*alp(i)*dsqrt((v1-v0)**2*sv1
     &         +(v2-v0)**2*sv2+(v1-v2)**2*sv0)
        enddo
C
        write(6,*) ''
        write(6,*) ''
        write(6,*) ''
        write(6,*) 'chi2    = ',chi2
        write(6,*) 'Cx1     = ',val(1),sval(1)
        write(6,*) 'Cx2     = ',val(2),sval(2)
        write(6,*) 'Cx0     = ',val(3),sval(3)
        write(6,*) 'Hx1     = ',val(4),sval(4)
        write(6,*) 'Hx2     = ',val(5),sval(5)
        write(6,*) 'Hx0     = ',val(6),sval(6)
        write(6,*) ''
        write(6,*) 'Cx0     = ',val(3),sval(3)
        write(6,*) 'Cg0     = ',alp(1)*(val(1)-val(3))**2,
     &              alp(1)*(val(1)-val(3))**2*dsqrt((salp(1)/alp(1))**2
     &              +4.d0/(val(1)-val(3))**2*(sval(1)**2+sval(3)**2)) 
        write(6,*) 'Calp    = ',alp(1),salp(1)
        write(6,*) 'Hx0     = ',val(6),sval(6)
        write(6,*) 'Hg0     = ',alp(2)*(val(4)-val(6))**2,
     &              alp(2)*(val(4)-val(6))**2*dsqrt((salp(2)/alp(2))**2
     &              +4.d0/(val(6)-val(4))**2*(sval(4)**2+sval(6)**2)) 
        write(6,*) 'Halp    = ',alp(2),salp(2)
        if(val(2).gt.val(3)) then
          write(6,*) 'Cin     =    1'
        else
          write(6,*) 'Cin     =   -1'
        endif
        if(val(5).gt.val(6)) then
          write(6,*) 'Hin     =    1'
        else
          write(6,*) 'Hin     =   -1'
        endif
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine inv(mat,val,nmat)
        integer nmax
        parameter(nmax=100)
        double precision mat(nmax,nmax)
        double precision val(nmax)
        integer ipiv(nmax)
        integer i,j,k,irow,icol,nmat,imat
        integer indxc(nmax),indxr(nmax)
        double precision big,pivinv,dum
C
        do i=1,nmat
          ipiv(i)=0
        enddo
C
        do i=1,nmat
          big=0.d0
          do j=1,nmat
            if(ipiv(j).eq.0) then
              do k=1,nmat
                if(ipiv(k).eq.0) then
                  if(abs(mat(j,k)).ge.big) then
                    big=abs(mat(j,k))
                    irow=j
                    icol=k
                  endif
                endif
              enddo
            endif
          enddo
          ipiv(icol)=ipiv(icol)+1
C
          if(irow.ne.icol) then
            do j=1,nmat
              dum=mat(irow,j)
              mat(irow,j)=mat(icol,j)
              mat(icol,j)=dum
            enddo
            dum=val(irow)
            val(irow)=val(icol)
            val(icol)=dum
          endif
          indxr(i)=irow
          indxc(i)=icol
          pivinv=1.d0/mat(icol,icol)
          mat(icol,icol)=1.d0
          do j=1,nmat
            mat(icol,j)=mat(icol,j)*pivinv
          enddo
          val(icol)=val(icol)*pivinv
C
          do j=1,nmat
            if(j.ne.icol) then
              dum=mat(j,icol)
              mat(j,icol)=0.d0
              do k=1,nmat
                mat(j,k)=mat(j,k)-mat(icol,k)*dum
              enddo
              val(j)=val(j)-val(icol)*dum
            endif
          enddo
        enddo
C
        do i=nmat,1,-1
          if(indxr(i).ne.indxc(i)) then
            do j=1,nmat
              dum=mat(j,indxr(i))
              mat(j,indxr(i))=mat(j,indxc(i))
              mat(j,indxc(i))=dum
            enddo
          endif
        enddo
C
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
