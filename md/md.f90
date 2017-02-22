



subroutine run_simulation
	use config
	use atom_prop
	implicit none
	integer step
	integer outsideSteps

	outsideSteps = nSteps/deltaWrite	! assuming they divide nicely

	call open_output_files()
	do step = 1, outsideSteps
		! open acc protocol? send all data
		call perform_md(deltaWrite)
		! end acc protocol? receive all data
		call write_output_files()

	enddo
	call close_output_files()
	

endsubroutine run_simulation


subroutine perform_md(stepsToPerform)
	use config
	use forcefield
	use atom_params
	use random
	use initialize
	implicit none
	integer stepsToPerform
	integer step
       

	do step=1, stepsToPerform 
!  Get new forces
		call getfrc()
!  Increment velocity and position
!$OMP parallel do schedule(static) num_threads(ncpu) private(i,j,try) shared(natom,x,v,f,pnu,T,mass,hdt,dt)
		do i=1,natom
			try=rand()
			if (try.lt.pnu) then
				call thermo(v(:,i),mass(i))
				do j=1,3
			    		v(j,i)=v(j,i)+f(j,i)/mass(i)*hdt
			    		x(j,i)=x(j,i)+dt*v(j,i)
			  	enddo
			else
				do j=1,3
					v(j,i)=v(j,i)+f(j,i)/mass(i)*dt
					x(j,i)=x(j,i)+dt*v(j,i)
				enddo
			endif
		enddo
!$OMP end parallel do
		call wrap()
	enddo
            call getcom(r2,rcom)
            write(15,999) igo,jgo-1,r2
999         format(i7,i1,1X,f16.12)

endsubroutine perform_md

subroutine open_output_files
	use rst_file
	use config, only : fname
	use atom_params
	implicit none
	character*80 xyz
	character*80 vel
	character*80 frc
	character*80 ene
	character*80 com

        write(rst,555) trim(fname)
555     format(a,'.rst')
        write(xyz,554) trim(fname)
554     format(a,'.mdcrd')
        write(vel,553) trim(fname)
553     format(a,'.mdvel')
        write(frc,552) trim(fname)
552     format(a,'.mdfrc')
        write(ene,551) trim(fname)
551     format(a,'.log')
        write(com,550) trim(fname)
550     format(a,'.dat')
C
        open(19,FILE=xyz,STATUS='unknown')
        write(19,'(a12)') "default_name"
        open(18,FILE=vel,STATUS='unknown')
        write(18,'(a12)') "default_name"
        open(17,FILE=frc,STATUS='unknown')
        write(17,'(a12)') "default_name"
        open(16,FILE=ene,STATUS='unknown')
        open(15,FILE=com,STATUS='unknown')
        call getene(x,v,lbox,hlbox,16,0)

endsubroutine open_output_files

subroutine write_output_files()
	use forcefield, only : natom
	use atom_params
	implicit none

          call writecrd(x,lbox,natom,19)
          call writecrd(v,lbox,natom,18)
          call writecrd(f,lbox,natom,17)

endsubroutine write_output_files()

subroutine close_output_files()
	use rst_file
	use atom_params
	use forcefield, only : natom
	implicit none

        call writerst(x,v,lbox,natom,rst)
        close(19)
        close(18)
        close(17)
        close(16)
        close(15)

endsubroutine close_output_files



subroutine getene(nf,nstep)
        use forcefield
	implicit none
        integer i,j,k,iex,jex,it,jt,nlj
        integer i1,i2,i3,i4,i5,nf,nstep
        double precision r2,r6,elj,ec,ebnd,eang,e14e,e14v,edih,ekin,erst
        double precision r12,r23,rdot,th,d1(3),d2(3),rcom(3)
        double precision c11,c12,c13,c22,c23,c33,a,b

        call getcom(x,lbox,hlbox,umblist,mass,r2,rcom)
        erst=10.d0*(r2-rumb)**2

        elj=0.d0
        ec=0.d0
        iex=1
        if (excl(1).eq.0) then
          jex=1
        else
          jex=0
        endif
        do i=1,natom
          it=ityp(i)
          do j=i+1,natom
            if (excl(iex).ne.j) then
              jt=ityp(j)
              nlj=ntyp*(it-1)+jt
              nlj=nbparm(nlj)
              r2=0.d0
              do k=1,3
                dx=x(k,i)-x(k,j)
                if (dx.gt.hlbox) then
                  dx=dx-lbox
                endif
                if (dx.lt.-hlbox) then
                  dx=dx+lbox
                endif
                r2=r2+dx*dx
              enddo
              r6=r2**(-3)
              elj=elj+r6*(alj(nlj)*r6-blj(nlj))
              ec=ec+charge(i)*charge(j)/dsqrt(r2)
            else
              iex=iex+1
            endif
          enddo
          if (excl(iex).eq.0) then
            if (jex.eq.0) then
              jex=1
            else
              iex=iex+1
              if (excl(iex).eq.0) then
                jex=1
              else
                jex=0
              endif
            endif
          endif
        enddo

        ebnd=0.d0
        do i=1,nbona
          i1=balist(1,i)
          i2=balist(2,i)
          it=balist(3,i)
          r2=0.d0
          do k=1,3
            r2=r2+(x(k,i1)-x(k,i2))**2
          enddo
          ebnd=ebnd+kbnd(it)*(dsqrt(r2)-rbnd(it))**2
        enddo
        do i=1,nbonh
          i1=bhlist(1,i)
          i2=bhlist(2,i)
          it=bhlist(3,i)
          r2=0.d0
          do k=1,3
            r2=r2+(x(k,i1)-x(k,i2))**2
          enddo
          ebnd=ebnd+kbnd(it)*(dsqrt(r2)-rbnd(it))**2
        enddo

        eang=0.d0
        do i=1,ntheta
          i1=aalist(1,i)
          i2=aalist(2,i)
          i3=aalist(3,i)
          it=aalist(4,i)
          c11=0.d0
          c22=0.d0
          c12=0.d0
          do k=1,3
            d1(k)=x(k,i1)-x(k,i2)
            d2(k)=x(k,i2)-x(k,i3)
            c11=c11+d1(k)*d1(k)
            c22=c22+d2(k)*d2(k)
            c12=c12+d1(k)*d2(k)
          enddo
          b=-c12/dsqrt(c11*c22)
          if (b.ge.1.d0) then
            th=1.d-16
          elseif (b.le.-1.d0) then
            th=3.1415926535d0
          else
            th=dacos(b)
          endif
          eang=eang+kang(it)*(th-tang(it))**2
        enddo
        do i=1,ntheth
          i1=ahlist(1,i)
          i2=ahlist(2,i)
          i3=ahlist(3,i)
          it=ahlist(4,i)
          c11=0.d0
          c22=0.d0
          c12=0.d0
          do k=1,3
            d1(k)=x(k,i1)-x(k,i2)
            d2(k)=x(k,i2)-x(k,i3)
            c11=c11+d1(k)*d1(k)
            c22=c22+d2(k)*d2(k)
            c12=c12+d1(k)*d2(k)
          enddo
          b=-c12/dsqrt(c11*c22)
          if (b.ge.1.d0) then
            th=1.d-16
          elseif (b.le.-1.d0) then
            th=3.1415926535d0
          else
            th=dacos(b)
          endif
          eang=eang+kang(it)*(th-tang(it))**2
        enddo

        edih=0.d0
        e14e=0.d0
        e14v=0.d0
        do i=1,nphih
          i1=dhlist(1,i)
          i2=dhlist(2,i)
          i3=dhlist(3,i)
          i4=dhlist(4,i)
          i5=dhlist(5,i)
          if (i4.lt.0) i4=-i4
          if (i3.lt.0) then
            i3=-i3
          else
            r2=0.d0
            do k=1,3
              r2=r2+(x(k,i1)-x(k,i4))*(x(k,i1)-x(k,i4))
            enddo
            r6=r2**(-3)
            it=ityp(i1)
            jt=ityp(i4)
            nlj=ntyp*(it-1)+jt
            nlj=nbparm(nlj)
            e14e=e14e+charge(i1)*charge(i4)/dsqrt(r2)/scee(i5)
            e14v=e14v+r6*(alj(nlj)*r6-blj(nlj))/scnb(i5)
          endif

          if (i2.gt.0) then
            c11=0.d0
            c12=0.d0
            c13=0.d0
            c22=0.d0
            c23=0.d0
            c33=0.d0
            do k=1,3
              c11=c11+(x(k,i1)-x(k,i2))*(x(k,i1)-x(k,i2))
              c12=c12+(x(k,i1)-x(k,i2))*(x(k,i2)-x(k,i3))
              c13=c13+(x(k,i1)-x(k,i2))*(x(k,i3)-x(k,i4))
              c22=c22+(x(k,i2)-x(k,i3))*(x(k,i2)-x(k,i3))
              c23=c23+(x(k,i2)-x(k,i3))*(x(k,i3)-x(k,i4))
              c33=c33+(x(k,i3)-x(k,i4))*(x(k,i3)-x(k,i4))
            enddo
            a=c12*c23-c13*c22
            b=(c11*c22-c12*c12)*(c22*c33-c23*c23)
            b=a/dsqrt(b)
            if (b.le.-1.d0) then
              th=3.1415926535d0
            elseif (b.ge.1.d0) then
              th=1.d-16
            else
              th=dacos(b)
            endif
            edih=edih+kdih(i5)*(1.d0+dcos(ndih(i5)*th-pdih(i5)))
          else
            edih=edih+kdih(i5)*(1.d0+dcos(ndih(i5)*th-pdih(i5)))
          endif
        enddo
        do i=1,nphia
          i1=dalist(1,i)
          i2=dalist(2,i)
          i3=dalist(3,i)
          i4=dalist(4,i)
          i5=dalist(5,i)
          if (i4.lt.0) i4=-i4
          if (i3.lt.0) then
            i3=-i3
          else
            r2=0.d0
            do k=1,3
              r2=r2+(x(k,i1)-x(k,i4))*(x(k,i1)-x(k,i4))
            enddo
            r6=r2**(-3)
            it=ityp(i1)
            jt=ityp(i4)
            nlj=ntyp*(it-1)+jt
            nlj=nbparm(nlj)
            e14e=e14e+charge(i1)*charge(i4)/dsqrt(r2)/scee(i5)
            e14v=e14v+r6*(alj(nlj)*r6-blj(nlj))/scnb(i5)
          endif

          if (i2.gt.0) then
            c11=0.d0
            c12=0.d0
            c13=0.d0
            c22=0.d0
            c23=0.d0
            c33=0.d0
            do k=1,3
              c11=c11+(x(k,i1)-x(k,i2))*(x(k,i1)-x(k,i2))
              c12=c12+(x(k,i1)-x(k,i2))*(x(k,i2)-x(k,i3))
              c13=c13+(x(k,i1)-x(k,i2))*(x(k,i3)-x(k,i4))
              c22=c22+(x(k,i2)-x(k,i3))*(x(k,i2)-x(k,i3))
              c23=c23+(x(k,i2)-x(k,i3))*(x(k,i3)-x(k,i4))
              c33=c33+(x(k,i3)-x(k,i4))*(x(k,i3)-x(k,i4))
            enddo
            a=c12*c23-c13*c22
            b=(c11*c22-c12*c12)*(c22*c33-c23*c23)
            b=a/dsqrt(b)
            if (b.le.-1.d0) then
              th=3.1415926535d0
            elseif (b.ge.1.d0) then
              th=1.d-16
            else
              th=dacos(b)
            endif
            edih=edih+kdih(i5)*(1.d0+dcos(ndih(i5)*th-pdih(i5)))
          else
            edih=edih+kdih(i5)*(1.d0+dcos(ndih(i5)*th-pdih(i5)))
          endif
        enddo

        ekin=0.d0
        do i=1,natom
          do k=1,3
            ekin=ekin+0.5d0*mass(i)*v(k,i)*v(k,i)
          enddo
        enddo

        write(nf,997) nstep
        write(nf,999) "ETOT      =   ",elj+ec+ebnd+eang+edih+erst &
             +e14v+e14e+ekin,"KINETIC   =   ",ekin,"POTENTIAL =   ", &
             elj+ec+ebnd+eang+edih+e14v+e14e+erst
        write(nf,999) "BOND      =   ",ebnd,"ANGLE     =   ",eang, &
             "DIHEDRAL  =   ",edih
        write(nf,998) "VDWAALS   =   ",elj,"EEL       =   ",ec
        write(nf,998) "1-4 VDW   =   ",e14v,"1-4 EEL   =   ",e14e
        write(nf,996) "RESTRAINT =   ",erst
        write(nf,*) ''
999     format(3(2X,a14,f8.4))
998     format(2(2X,a14,f8.4))
996     format(2X,a14,f8.4)
997     format('TIME STEP:  ',i8)
C
      return
      end

subroutine getfrc()
        use atom_params
	use forcefield
	use random
!	use solvent
	implicit none

        integer i,j,k,iex,jex,it,jt,nlj
        integer i1,i2,i3,i4,i5
        double precision flj,fc,fbnd,fang,f14e,f14v,fdih,frst1,frst2
        double precision dis(3,natom)
        double precision r2,r6,r12,r23,rdot,th,r(3),d1(3),d2(3),d3(3)
        double precision c11,c12,c13,c22,c23,c33,a,b
        double precision t1,t2,t3,t4,t5,t6,f1,f4
!
        double precision rnow,x1,x2,x3,y1,y2,y3,dx,dy,dz
        double precision gnow,fs,try,prob
        double precision rand
        double precision rcom(3)
        integer igo,ngo

!!$OMP threadprivate(/random/)
!!
!!  Solvent!
!!$OMP parallel do schedule(dynamic) num_threads(ncpu)
!!$OMP& private(i,it,j,jt,k,ngo,rnow,prob,try,x1,x2,x3,y1,y2,y3,r2
!!$OMP&             ,gnow,igo,dx,dy,dz,r6,fs,dis)
!!$OMP& shared(natom,npts,w,x0,ityp,gr2,alp,g0,f,x,lbox,hlbox)
	do i=1,natom
		do j=1,3
			f(j,i)=0.d0
		enddo
		it=ityp(i)
		do j=1,natom
			if (j.ne.i) then
				do k=1,3
					dis(k,j)=x(k,i)-x(k,j)
					if (dis(k,j).gt.hlbox) then
						dis(k,j)=dis(k,j)-lbox
					endif
					if (dis(k,j).lt.-hlbox) then
						dis(k,j)=dis(k,j)+lbox
					endif
				enddo
			endif
		enddo
		do ngo=1,npts
			rnow=1.d0-2.d0*rand()
			prob=rnow*rnow
			try=rand()
			do while(try.lt.prob)
				rnow=1.d0-2.d0*rand()
				prob=rnow*rnow
				try=rand()
			enddo
			rnow=w(it)*rnow+x0(it)
			x1=1.d0-2.d0*rand()
			x2=1.d0-2.d0*rand()
			r2=x1*x1+x2*x2
			do while (r2.gt.1.d0)
				x1=1.d0-2.d0*rand()
				x2=1.d0-2.d0*rand()
				r2=x1*x1+x2*x2
			enddo
			y1=rnow*(1.d0-2.d0*r2)
			r2=2.d0*dsqrt(1.d0-r2)
			y2=rnow*x1*r2
			y3=rnow*x2*r2
			
			gnow=1.d0
			igo=1
			j=0
			do while(igo.eq.1)
				j=j+1
				if (j.eq.i) j=j+1
				if (i.eq.natom) then
					if (j.eq.natom-1) igo=-1
				else
					if (j.eq.natom) igo=-1
				endif
				jt=ityp(j)
				
				dx=y1+dis(1,j)
				dy=y2+dis(2,j)
				dz=y3+dis(3,j)
				r2=dx*dx+dy*dy+dz*dz
				if (r2.lt.gr2(1,jt)) then
					igo=0
				elseif (r2.lt.gr2(2,jt)) then
					gnow=gnow*(-alp(jt)*(dsqrt(r2)-x0(jt))**2+g0(jt))
				endif
			enddo
			
			if (igo.eq.-1) then
				r6=rnow**(-6)
				fs=gnow*r6*(Bs(it)-As(it)*r6)
				f(1,i)=f(1,i)+fs*y1
				f(2,i)=f(2,i)+fs*y2
				f(3,i)=f(3,i)+fs*y3
			endif
		enddo
		do j=1,3
			f(j,i)=f(j,i)*vtot(it)
		enddo
	enddo
!!C$OMP end parallel do
!  Restraint
        call getcom(x,lbox,hlbox,umblist,mass,r2,rcom)
        frst1=kumb1*(r2-rumb)
        frst2=kumb2*(r2-rumb)
        do i=2,umblist(1,1)
          i1=umblist(i,1)
          do j=1,3
            f(j,i1)=f(j,i1)-frst1*rcom(j)*mass(i1)
          enddo
        enddo
        do i=2,umblist(1,2)
          i1=umblist(i,2)
          do j=1,3
            f(j,i1)=f(j,i1)+frst2*rcom(j)*mass(i1)
          enddo
        enddo
!  Non-bonded
        iex=1
        if (excl(1).eq.0) then
          jex=1
        else
          jex=0
        endif
        do i=1,natom
          it=ityp(i)
          do j=i+1,natom
            if (excl(iex).ne.j) then
              jt=ityp(j)
              nlj=ntyp*(it-1)+jt
              nlj=nbparm(nlj)
              r2=0.d0
              do k=1,3
                dx=x(k,i)-x(k,j)
                if (dx.gt.hlbox) then
                  dx=dx-lbox
                endif
                if (dx.lt.-hlbox) then
                  dx=dx+lbox
                endif
                r(k)=dx
                r2=r2+dx*dx
              enddo
              r6=r2**(-3)
              flj=r6*(12.d0*alj(nlj)*r6-6.d0*blj(nlj))/r2
              fc=charge(i)*charge(j)/dsqrt(r2)/r2
              do k=1,3
                f(k,i)=f(k,i)+(flj+fc)*r(k)
                f(k,j)=f(k,j)-(flj+fc)*r(k)
              enddo
            else
              iex=iex+1
            endif
          enddo
          if (excl(iex).eq.0) then
            if (jex.eq.0) then
              jex=1
            else
              iex=iex+1
              if (excl(iex).eq.0) then
                jex=1
              else
                jex=0
              endif
            endif
          endif
        enddo
!  Bonds
        do i=1,nbona
          i1=balist(1,i)
          i2=balist(2,i)
          it=balist(3,i)
          r2=0.d0
          do k=1,3
            r(k)=x(k,i1)-x(k,i2)
            r2=r2+r(k)*r(k)
          enddo
          fbnd=2.d0*kbnd(it)*(rbnd(it)/dsqrt(r2)-1.d0)
          do k=1,3
            f(k,i1)=f(k,i1)+fbnd*r(k)
            f(k,i2)=f(k,i2)-fbnd*r(k)
          enddo
        enddo
        do i=1,nbonh
          i1=bhlist(1,i)
          i2=bhlist(2,i)
          it=bhlist(3,i)
          r2=0.d0
          do k=1,3
            r(k)=x(k,i1)-x(k,i2)
            r2=r2+r(k)**2
          enddo
          fbnd=2.d0*kbnd(it)*(rbnd(it)/dsqrt(r2)-1.d0)
          do k=1,3
            f(k,i1)=f(k,i1)+fbnd*r(k)
            f(k,i2)=f(k,i2)-fbnd*r(k)
          enddo
        enddo
!  Angles
        do i=1,ntheta
          i1=aalist(1,i)
          i2=aalist(2,i)
          i3=aalist(3,i)
          it=aalist(4,i)
          c11=0.d0
          c22=0.d0
          c12=0.d0
          do k=1,3
            d1(k)=x(k,i1)-x(k,i2)
            d2(k)=x(k,i2)-x(k,i3)
            c11=c11+d1(k)*d1(k)
            c22=c22+d2(k)*d2(k)
            c12=c12+d1(k)*d2(k)
          enddo
          b=-c12/dsqrt(c11*c22)
          if (b.ge.1.d0) then
            th=1.d-16
          elseif (b.le.-1.d0) then
            th=3.1415926535d0
          else
            th=dacos(b)
          endif
          fang=2.d0*kang(it)*(th-tang(it))/dsqrt(c11*c22-c12*c12)
          do k=1,3
            f(k,i1)=f(k,i1)+fang*(c12/c11*d1(k)-d2(k))
            f(k,i2)=f(k,i2)+fang*((1.d0+c12/c22)*d2(k)-(1.d0+c12/c11)*d1(k))
            f(k,i3)=f(k,i3)+fang*(d1(k)-c12/c22*d2(k))
          enddo
        enddo
        do i=1,ntheth
          i1=ahlist(1,i)
          i2=ahlist(2,i)
          i3=ahlist(3,i)
          it=ahlist(4,i)
          c11=0.d0
          c22=0.d0
          c12=0.d0
          do k=1,3
            d1(k)=x(k,i1)-x(k,i2)
            d2(k)=x(k,i2)-x(k,i3)
            c11=c11+d1(k)*d1(k)
            c22=c22+d2(k)*d2(k)
            c12=c12+d1(k)*d2(k)
          enddo
          b=-c12/dsqrt(c11*c22)
          if (b.ge.1.d0) then
            th=1.d-16
          elseif (b.le.-1.d0) then
            th=3.1415926535d0
          else
            th=dacos(b)
          endif
          fang=2.d0*kang(it)*(th-tang(it))/dsqrt(c11*c22-c12*c12)
          do k=1,3
            f(k,i1)=f(k,i1)+fang*(c12/c11*d1(k)-d2(k))
            f(k,i2)=f(k,i2)+fang*((1.d0+c12/c22)*d2(k)-(1.d0+c12/c11)*d1(k))
            f(k,i3)=f(k,i3)+fang*(d1(k)-c12/c22*d2(k))
          enddo
        enddo
!! Dihedrals
        do i=1,nphih
          i1=dhlist(1,i)
          i2=dhlist(2,i)
          i3=dhlist(3,i)
          i4=dhlist(4,i)
          i5=dhlist(5,i)
          if (i4.lt.0) i4=-i4
          if (i3.lt.0) then
            i3=-i3
          else
            r2=0.d0
            do k=1,3
              r(k)=x(k,i1)-x(k,i4)
              r2=r2+r(k)*r(k)
            enddo
            r6=r2**(-3)
            it=ityp(i1)
            jt=ityp(i4)
            nlj=ntyp*(it-1)+jt
            nlj=nbparm(nlj)
            f14e=charge(i1)*charge(i4)/r2/dsqrt(r2)/scee(i5)
            f14v=r6*(12.d0*alj(nlj)*r6-6.d0*blj(nlj))/scnb(i5)/r2
            do k=1,3
              f(k,i1)=f(k,i1)+(f14e+f14v)*r(k)
              f(k,i4)=f(k,i4)-(f14e+f14v)*r(k)
            enddo
          endif
!!!
          if (i2.gt.0) then
            c11=0.d0
            c12=0.d0
            c13=0.d0
            c22=0.d0
            c23=0.d0
            c33=0.d0
            do k=1,3
              d1(k)=x(k,i1)-x(k,i2)
              d2(k)=x(k,i2)-x(k,i3)
              d3(k)=x(k,i3)-x(k,i4)
              c11=c11+d1(k)*d1(k)
              c12=c12+d1(k)*d2(k)
              c13=c13+d1(k)*d3(k)
              c22=c22+d2(k)*d2(k)
              c23=c23+d2(k)*d3(k)
              c33=c33+d3(k)*d3(k)
            enddo
!
            t1=c13*c22-c12*c23
            t2=c11*c23-c12*c13
            t3=c12*c12-c11*c22
            t4=c22*c33-c23*c23
            t5=c13*c23-c12*c33
            t6=-t1
!
            b=dsqrt(-t3*t4)
!
            a=t6/b
            if (a.le.-1.d0) then
              fdih=0.d0
            elseif (a.ge.1.d0) then
              fdih=0.d0
            else
              th=dacos(a)
              fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))/dsin(th)*c22/b
            endif
          else
            i2=-i2
            if (a.le.-1.d0) then
              fdih=0.d0
            elseif (a.ge.1.d0) then
              fdih=0.d0
            else
              th=dacos(a)
              fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))/dsin(th)*c22/b
            endif
          endif
          do k=1,3
            f1=fdih*(t1*d1(k)+t2*d2(k)+t3*d3(k))/t3
            f4=-fdih*(t4*d1(k)+t5*d2(k)+t6*d3(k))/t4
            f(k,i1)=f(k,i1)+f1
            f(k,i2)=f(k,i2)-(1.d0+c12/c22)*f1+c23/c22*f4
            f(k,i3)=f(k,i3)+c12/c22*f1-(1.d0+c23/c22)*f4
            f(k,i4)=f(k,i4)+f4
          enddo
        enddo
        do i=1,nphia
          i1=dalist(1,i)
          i2=dalist(2,i)
          i3=dalist(3,i)
          i4=dalist(4,i)
          i5=dalist(5,i)
          if (i4.lt.0) i4=-i4
          if (i3.lt.0) then
            i3=-i3
          else
            r2=0.d0
            do k=1,3
              r(k)=x(k,i1)-x(k,i4)
              r2=r2+r(k)*r(k)
            enddo
            r6=r2**(-3)
            it=ityp(i1)
            jt=ityp(i4)
            nlj=ntyp*(it-1)+jt
            nlj=nbparm(nlj)
            f14e=charge(i1)*charge(i4)/r2/dsqrt(r2)/scee(i5)
            f14v=r6*(12.d0*alj(nlj)*r6-6.d0*blj(nlj))/scnb(i5)/r2
            do k=1,3
              f(k,i1)=f(k,i1)+(f14e+f14v)*r(k)
              f(k,i4)=f(k,i4)-(f14e+f14v)*r(k)
            enddo
          endif
!
          if (i2.gt.0) then
            c11=0.d0
            c12=0.d0
            c13=0.d0
            c22=0.d0
            c23=0.d0
            c33=0.d0
            do k=1,3
              d1(k)=x(k,i1)-x(k,i2)
              d2(k)=x(k,i2)-x(k,i3)
              d3(k)=x(k,i3)-x(k,i4)
              c11=c11+d1(k)*d1(k)
              c12=c12+d1(k)*d2(k)
              c13=c13+d1(k)*d3(k)
              c22=c22+d2(k)*d2(k)
              c23=c23+d2(k)*d3(k)
              c33=c33+d3(k)*d3(k)
            enddo
!
            t1=c13*c22-c12*c23
            t2=c11*c23-c12*c13
            t3=c12*c12-c11*c22
            t4=c22*c33-c23*c23
            t5=c13*c23-c12*c33
            t6=-t1
!
            b=dsqrt(-t3*t4)
!
            a=t6/b
            if (a.le.-1.d0) then
              fdih=0.d0
            elseif (a.ge.1.d0) then
              fdih=0.d0
            else
              th=dacos(a)
              fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))/dsin(th)*c22/b
            endif
          else
            i2=-i2
            if (a.le.-1.d0) then
              fdih=0.d0
            elseif (a.ge.1.d0) then
              fdih=0.d0
            else
              th=dacos(a)
              fdih=ndih(i5)*kdih(i5)*dsin(ndih(i5)*th-pdih(i5))/dsin(th)*c22/b
            endif
          endif
          do k=1,3
            f1=fdih*(t1*d1(k)+t2*d2(k)+t3*d3(k))/t3
            f4=-fdih*(t4*d1(k)+t5*d2(k)+t6*d3(k))/t4
            f(k,i1)=f(k,i1)+f1
            f(k,i2)=f(k,i2)-(1.d0+c12/c22)*f1+c23/c22*f4
            f(k,i3)=f(k,i3)+c12/c22*f1-(1.d0+c23/c22)*f4
            f(k,i4)=f(k,i4)+f4
          enddo
        enddo
!
      return
      end
!
subroutine writerst(x,v,lbox,natom,rst)
        integer nmax
        parameter(nmax=10000)
        integer natom
        double precision x(3,nmax),v(3,nmax)
        double precision lbox
        character*64 rst
        integer i,ip
C
        open(20,FILE=rst,STATUS='unknown')
        write(20,999)
        write(20,998) natom,1.d0
        i=1
        do ip=1,natom/2
          write(20,997) x(1,i),x(2,i),x(3,i),x(1,i+1),x(2,i+1),x(3,i+1)
          i=i+2
        enddo
        if (i-1.ne.natom) then
          write(20,996) x(1,natom),x(2,natom),x(3,natom)
        endif
        i=1
        do ip=1,natom/2
          write(20,997) v(1,i),v(2,i),v(3,i),v(1,i+1),v(2,i+1),v(3,i+1)
          i=i+2
        enddo
        if (i-1.ne.natom) then
          write(20,996) v(1,natom),v(2,natom),v(3,natom)
        endif
        write(20,997) lbox,lbox,lbox,90.d0,90.d0,90.d0
        close(20)
!
999     format('Rex Generated Restart')
998     format(i5,e15.7)
997     format(6(f12.7))
996     format(6(f12.7))
!
      return
      end
!
subroutine writecrd(x,lbox,natom,nf)
	implicit none
        integer natom
        double precision x(3,natom)
        double precision lbox
        integer nf
!
        integer i,j,k
!
        k=1
        do i=1,natom
          do j=1,3
            write(nf,999) x(j,i)
            k=k+1
            if (k.eq.11) then
              write(nf,998)
              k=1
            endif
          enddo
        enddo
        if (k.ne.11) write(nf,998)
        if (nf.eq.19) then
          do i=1,3
            write(nf,999) lbox
          enddo
          write(nf,998)
        endif
999     format(f8.3,$)
998     format('')
!
      return
endsubroutine writecrd
!
subroutine getcom(r2,rcom)
	use forcefield, only : natom, mass
	use umbrella
	use atom_params
	implicit none
        double precision xnow(3)
        double precision xcom(3,2),mtot,r2,dx
        double precision rcom(3)
!
        integer i,j,i1,i2,nf,ip,i0
!
        do j=1,3
          xcom(j,1)=0.d0
          xcom(j,2)=0.d0
        enddo
        mtot=0.d0
        i0=umSel1(2)
        do i=2,umSel1(1)
          ip=umSel1(i)
          do j=1,3
            dx=x(j,ip)-x(j,i0)
            if (dx.gt.hlbox) then
              xnow(j)=x(j,ip)-lbox
            elseif (dx.lt.-hlbox) then
              xnow(j)=x(j,ip)+lbox
            else
              xnow(j)=x(j,ip)
            endif
          enddo
          do j=1,3
            xcom(j,1)=xcom(j,1)+mass(ip)*xnow(j)
          enddo
          mtot=mtot+mass(ip)
        enddo
        do j=1,3
          xcom(j,1)=xcom(j,1)/mtot
        enddo
C
        mtot=0.d0
        i0=umblist(2,2)
        do i=2,umblist(1,2)
          ip=umblist(i,2)
          do j=1,3
            dx=x(j,ip)-x(j,i0)
            if (dx.gt.hlbox) then
              xnow(j)=x(j,ip)-lbox
            elseif (dx.lt.-hlbox) then
              xnow(j)=x(j,ip)+lbox
            else
              xnow(j)=x(j,ip)
            endif
          enddo
          do j=1,3
            xcom(j,2)=xcom(j,2)+mass(ip)*xnow(j)
          enddo
          mtot=mtot+mass(ip)
        enddo
        do j=1,3
          xcom(j,2)=xcom(j,2)/mtot
        enddo
C
        r2=0.d0
        do i=1,3
          dx=xcom(i,1)-xcom(i,2)
          if (dx.gt.hlbox) then
            dx=dx-lbox
          endif
          if (dx.lt.-hlbox) then
            dx=dx+lbox
          endif
          rcom(i)=dx
          r2=r2+dx*dx
        enddo
        r2=dsqrt(r2)
        do i=1,3
          rcom(i)=rcom(i)/r2
        enddo
C
      return
endsubroutine getcom
