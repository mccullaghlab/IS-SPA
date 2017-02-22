
module initialize


	contains

	! read all input data to initialize simulation
	subroutine initialize_run
		use forcefield
		use config
		use omp
        	character*64 frst, prmtop, fname
		! read input
		call getinput(prmtop,frst,fname)
		! set parameters
		ncpu=2
		dt=0.002d0
		T=298.d0
		
		dt=dt*20.455d0
		hdt=dt*0.5d0
		T=T*0.00198717d0
		pnu=0.001d0
		! read prmtop file
        	call setparam(prmtop)
        	call makemol()
		!  Get the positions, velocities, and atom types
        	call getx(frst)
		! wrap coordinates
        	call wrap()

	endsubroutine initialize_run


	! read the input information from a text file
        subroutine getinput(prmtop,frst,fname)
		use config
		use umbrella
		implicit none
        	character*64 frst,fname,prmtop
        	integer numb1,numb2,i
!
		open(20,FILE="input.dat",STATUS='old')
		read(20,999) prmtop   ! amber formatted prmtop file
		read(20,999) frst     ! restart file name
		read(20,*) ivel       ! initialize velocities
		read(20,999) fname    ! ?
		read(20,*) rumb       ! umbrella sampling distance value
		read(20,*) numb1      ! number of atoms in selection 1
		allocate(umSel1(numb1+1))
		umSel1(1)=numb1+1
		do i=1,numb1
		  read(20,*) umSel1(i+1)
		enddo
		read(20,*) numb2
		allocate(umSel2(numb2+1))
		umSel2(1)=numb2+1
		do i=1,numb2
		  read(20,*) umSel2(i+1)
		enddo
999     format(a)
        	close(20)
!
      		return
      endsubroutine getinput

     ! read parmtop file
      subroutine setparam(prmtop)
	use forcefield
	implicit none
        integer natyp
        integer i,j,inext,i1,i2,i3,i4
        character*80 str, temp
        character*64 prmtop
	integer counter
!
!
        open(22,FILE=prmtop,STATUS='old')
!  numbers
        do i=1,6
          read(22,*)
        enddo
        read(22,999) natom,ntyp,nbonh,nbona,ntheth,ntheta,nphih,nphia
        read(22,998) nnb,nres,numbnd,numang,nptra,natyp
!  charge
	allocate(charge(natom))
        inext=7+(natom-1)/20
        do i=1,inext
          read(22,*)
        enddo
        inext=(natom-1)/5+1
	do i = 1, natom
		read(22,advance='no', '(e16.8)') charge(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  mass
	allocate(mass(natom))
        inext=5+int(float(natom-1)/10.0)
        do i=1,inext
          read(22,*)
        enddo
        inext=(natom-1)/5+1
	do i = 1, natom
		read(22,advance='no', '(e16.8)') mass(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  ityp ATOM_TYPE_INDEX
	allocate(ityp(natom))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, natom
		read(22,advance='no', '(i8)') ityp(i)
		if (mod(i,10)==0) read(22,*)
	enddo
!  n-excluded NUMBER_EXCLUDED_ATOMS
	allocate(nexc(natom))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, natom
		read(22,advance='no', '(i8)') nexc(i)
		if (mod(i,10)==0) read(22,*)
	enddo
! non-bonded parm NONBONDED_PARM_INDEX
	allocate(nbparm(ntyp*ntyp))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, ntyp*ntyp
		read(22,advance='no', '(i8)') nbparm(i)
		if (mod(i,10)==0) read(22,*)
	enddo
!  residue pointers RESIDUE_POINTER
	allocate(nrpnt(nres+1))
        inext=5+(nres-1)/20+1
        do i=1,inext
          read(22,*) 
        enddo
	do i = 1, nres
		read(22,advance='no', '(i8)') nrpnt(i)
		if (mod(i,10)==0) read(22,*)
	enddo
        nrpnt(nres+1)=natom+1
!  bond force
	allocate(kbnd(numbnd))
	inext=2
        do i=1,inext
          read(22,*) 
        enddo
	do i = 1, numbnd
		read(22,advance='no', '(e16.8)') kbnd(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  bond distance
	allocate(rbnd(numbnd))
        inext=3
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, numbnd
		read(22,advance='no', '(e16.8)') rbnd(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  angle force
	allocate(kang(numang))
        inext=3
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, numang
		read(22,advance='no', '(e16.8)') kang(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  angle values
	allocate(tang(numang))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, numang
		read(22,advance='no', '(e16.8)') tang(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  dihedral force
	allocate(kdih(nptra))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, nptra
		read(22,advance='no', '(e16.8)') kdih(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  dihedral preriodicity
	allocate(ndih(nptra))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, nptra
		read(22,advance='no', '(e16.8)') ndih(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  dihedral phase
	allocate(pdih(nptra))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, nptra
		read(22,advance='no', '(e16.8)') pdih(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  SCEE scale factor
	allocate(scee(nptra))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, nptra
		read(22,advance='no', '(e16.8)') scee(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  SCNB scale factor
	allocate(scnb(nptra))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, nptra
		read(22,advance='no', '(e16.8)') scnb(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  Lennard-Jones A
	allocate(alj(ntyp*(ntyp+1)/2))
        inext=5+(natyp-1)/5
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, ntyp*(ntyp+1)/2
		read(22,advance='no', '(e16.8)') alj(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  Lennard-Jones B
	allocate(blj(ntyp*(ntyp+1)/2))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i = 1, ntyp*(ntyp+1)/2
		read(22,advance='no', '(e16.8)') blj(i)
		if (mod(i,5)==0) read(22,*)
	enddo
!  Bond list w/ H
	allocate(bhlist(3,nbonh))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	counter = 0
	do i=1, nbonh
		do j=1,3
			counter = counter + 1
			read(22,advance='no','(i8)') bhlist(j,i)
			if (mod(counter,10)==0) read(22,*)
		enddo
	enddo
	! convert from coordinate array index to atom index
        do i=1,nbonh
          bhlist(1,i)=bhlist(1,i)/3+1
          bhlist(2,i)=bhlist(2,i)/3+1
        enddo
!  Bond list w/o H
	allocate(balist(3,nbona))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	counter = 0
	do i=1, nbona
		do j=1,3
			counter = counter + 1
			read(22,advance='no','(i8)') balist(j,i)
			if (mod(counter,10)==0) read(22,*)
		enddo
	enddo
	! convert from coordinate array index to atom index
        do i=1,nbona
          balist(1,i)=balist(1,i)/3+1
          balist(2,i)=balist(2,i)/3+1
        enddo
!  Angle list w/ H
	allocate(ahlist(4,ntheth))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	counter = 0
	do i=1, ntheth
		do j=1,4
			counter = counter + 1
			read(22,advance='no','(i8)') ahlist(j,i)
			if (mod(counter,10)==0) read(22,*)
		enddo
	enddo
	! convert from coordinate array index to atom index
        do i=1,ntheth
          ahlist(1,i)=ahlist(1,i)/3+1
          ahlist(2,i)=ahlist(2,i)/3+1
          ahlist(3,i)=ahlist(3,i)/3+1
        enddo
!  Angle list w/o H
	allocate(aalist(4,ntheta))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	counter = 0
	do i=1, ntheta
		do j=1,4
			counter = counter + 1
			read(22,advance='no','(i8)') aalist(j,i)
			if (mod(counter,10)==0) read(22,*)
		enddo
	enddo
	! convert from coordinate array index to atom index
        do i=1,ntheta
          aalist(1,i)=aalist(1,i)/3+1
          aalist(2,i)=aalist(2,i)/3+1
          aalist(3,i)=aalist(3,i)/3+1
        enddo
!  Dihedral list w/ H
	allocate(dhlist(5,nphih))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	counter = 0
	do i=1, nphih 
		do j=1,4
			counter = counter + 1
			read(22,advance='no','(i8)') dhlist(j,i)
			if (mod(counter,10)==0) read(22,*)
		enddo
	enddo
        do i=1,nphih
          dhlist(1,i)=dhlist(1,i)/3+1
          dhlist(2,i)=dhlist(2,i)/3+1
          if (dhlist(3,i).lt.0) then
            dhlist(3,i)=dhlist(3,i)/3-1
          else
            dhlist(3,i)=dhlist(3,i)/3+1
          endif
          if (dhlist(4,i).lt.0) then
            dhlist(4,i)=dhlist(4,i)/3-1
          else
            dhlist(4,i)=dhlist(4,i)/3+1
          endif
        enddo
        do i=2,nphih
          j=i-1
          if (iabs(dhlist(1,i)).eq.iabs(dhlist(1,j))) then
            if (iabs(dhlist(2,i)).eq.iabs(dhlist(2,j))) then
              if (iabs(dhlist(3,i)).eq.iabs(dhlist(3,j))) then
                if (iabs(dhlist(4,i)).eq.iabs(dhlist(4,j))) then
                  dhlist(2,i)=-dhlist(2,i)
                endif
              endif
            endif
          endif
        enddo
!  Dihedral list w/o H
	allocate(dalist(5,nphia))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	counter = 0
	do i=1, nphia 
		do j=1,4
			counter = counter + 1
			read(22,advance='no','(i8)') dalist(j,i)
			if (mod(counter,10)==0) read(22,*)
		enddo
	enddo
        do i=1,nphia
          dalist(1,i)=dalist(1,i)/3+1
          dalist(2,i)=dalist(2,i)/3+1
          if (dalist(3,i).lt.0) then
            dalist(3,i)=dalist(3,i)/3-1
          else
            dalist(3,i)=dalist(3,i)/3+1
          endif
          if (dalist(4,i).lt.0) then
            dalist(4,i)=dalist(4,i)/3-1
          else
            dalist(4,i)=dalist(4,i)/3+1
          endif
        enddo
        do i=2,nphia
          j=i-1
          if (iabs(dalist(1,i)).eq.iabs(dalist(1,j))) then
            if (iabs(dalist(2,i)).eq.iabs(dalist(2,j))) then
              if (iabs(dalist(3,i)).eq.iabs(dalist(3,j))) then
                if (iabs(dalist(4,i)).eq.iabs(dalist(4,j))) then
                  dalist(2,i)=-dalist(2,i)
                endif
              endif
            endif
          endif
        enddo
!  Excluded atom list
	allocate(excl(nnb))
        inext=2
        do i=1,inext
          read(22,*)
        enddo
	do i=1,nnb
		read(22,advance='no','(i8)') excl(i)
		if(mod(i,10)==0) read(22,*)
	enddo
!
999     format(8(i8))
998     format(2(i8),3(8X),4(i8))
997     format(5(e16.8))
996     format(10(i8))
        close(22)
!
      return
	endsubroutine

	! make neighbor list
      subroutine makemol()
	use forcefield
	use atom_prop
	implicit none!
        integer ires(natom),imol(nres),neigh(1000,1000)
        integer todo(nres)
        integer i,inow,i1,i2,j1,j2,n,k,j

	! allocate arrays
	allocate(x(3,natom))
	allocate(v(3,natom))
	allocate(f(3,natom))
!
        do i=1,nres
          neigh(1,i)=1
        enddo
!
        inow=0
        do i=1,natom
          if(i.eq.nrpnt(inow+1)) inow=inow+1
          ires(i)=inow
        enddo
        do i=1,nbonh
          i1=bhlist(1,i)
          i2=bhlist(2,i)
          j1=ires(i1)
          j2=ires(i2)
          if(j1.ne.j2) then
            n=neigh(1,j1)+1
            neigh(1,j1)=n
            neigh(n,j1)=j2
            n=neigh(1,j2)+1
            neigh(1,j2)=n
            neigh(n,j2)=j1
          endif
        enddo
        do i=1,nbona
          i1=balist(1,i)
          i2=balist(2,i)
          j1=ires(i1)
          j2=ires(i2)
          if(j1.ne.j2) then
            n=neigh(1,j1)+1
            neigh(1,j1)=n
            neigh(n,j1)=j2
            n=neigh(1,j2)+1
            neigh(1,j2)=n
            neigh(n,j2)=j1
          endif
        enddo
!
        nmol=0
        do i=1,nres
          imol(i)=0
        enddo
!
        do i=1,nres
          if (imol(i).eq.0) then
            nmol=nmol+1
            imol(i)=nmol
            i1=1
            i2=1
            todo(1)=i
            do while (i1.le.i2)
              j=todo(i1)
              j1=neigh(1,j)
              do k=2,j1
                j2=neigh(k,j)
                if (imol(j2).eq.0) then
                  i2=i2+1
                  todo(i2)=j2
                  imol(j2)=nmol
                endif
              enddo
              i1=i1+1
            enddo
          endif
        enddo
!
	allocate(nmpnt(nmol+1))
        j=0
        do i=1,nres
          if (imol(i).ne.j) then
            if (j+1.ne.imol(i)) then
              write(0,*) 'NONE SEQUENTIAL RESIDUES IN MOLECULE'
            endif
            j=j+1
            nmpnt(j)=nrpnt(i)
          endif
        enddo
        nmpnt(nmol+1)=natom+1
!
      return
      endsubroutine makemol
!!!!  
      subroutine getx(frst)
	use random
	use atom_prop
	use config, only : ivel
	use forcefield, only : mass, natom
	implicit none
	character*64 frst
        integer i,j,k,ip,ndx,jmin,jmax

        open(25,FILE=frst,STATUS='old')
        read(25,*)
        read(25,*)
        i=1
        ip=2
        do while (i.le.natom)
          read(25,999) x(1,i),x(2,i),x(3,i),x(1,ip),x(2,ip),x(3,ip)
          i=i+2
          ip=ip+2
        enddo
        if (ivel.eq.0) then
          do i=1,natom
            call thermo(v(:,i),mass(i))
          enddo
        else
          i=1
          ip=2
          do while (i.le.natom)
            read(25,999) v(1,i),v(2,i),v(3,i),v(1,ip),v(2,ip),v(3,ip)
            i=i+2
            ip=ip+2
          enddo
        endif
        read(25,998) lbox
999     format(6(f12.7))
998     format(f12.7)
        close(25)
        hlbox=lbox*0.5d0

      return
      endsubroutine getx

      subroutine wrap()
	use atom_prop
	use forcefield, only : nrpnt, nres
	implicit none
	double precision dx
        integer i,j,k,ndx,jmin,jmax

        do i=1,nres
          jmin=nrpnt(i)
          do k=1,3
            dx=x(k,jmin)
            if (dx.gt.lbox) then
              ndx=int(dx/lbox)
              jmax=nrpnt(i+1)-1
              do j=jmin,jmax
                x(k,j)=x(k,j)-dble(ndx)*lbox
              enddo
            endif
            if (dx.lt.0.d0) then
              ndx=int(-dx/lbox+1)
              jmax=nrpnt(i+1)-1
              do j=jmin,jmax
                x(k,j)=x(k,j)+dble(ndx)*lbox
              enddo
            endif
          enddo
        enddo
      return
      endsubroutine wrap

endmodule initialize

