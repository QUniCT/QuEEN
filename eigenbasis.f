        program eigenbasis
        implicit none
        integer::i,j,k,n,rn,info,ind
        double precision::time,psink,pground,normalize
        double precision,allocatable::a(:,:),w(:),work(:),rho(:)
        double complex,allocatable::r(:,:),re(:,:)
        double complex::im=(0,1)
        character(100)::frmt,fname,dm_name,h_name
        character(200),allocatable,dimension(:)::comlin
        logical::ex
!  open the files
        allocate(comlin(command_argument_count()))
        do i=1,command_argument_count(); call getarg(i,comlin(i)); enddo
        do i=1,command_argument_count()
         if (index(comlin(i),'-label').gt.0) read(comlin(i+1),'(a)') fname
        enddo
        write(dm_name,'(2a)')trim(fname),".dm"
        write(h_name,'(2a)')trim(fname),".hamiltonian"
        inquire(file=dm_name,exist=ex)
        if(.not.ex)stop"can't find the density matrix"
        inquire(file=h_name,exist=ex)
        if(.not.ex)stop"can't find the hamiltonian"
        open(60,file=h_name)
        open(61,file=dm_name)
!  calculate eigenvectors
        read(60,'(i10)') n
        allocate(a(n,n),w(n),work(3*n))
10      read(60,*,end=20) j,k,a(j,k); goto 10
20      do k=1,n; do j=1,k-1; a(k,j)=a(j,k); enddo; enddo
        call dsyev('V','U',n,a,n,w,work,3*n,info)
        write(frmt,'("(t5,i5,2x,es15.2,2x,",i0,"es10.2)")') n
!        do i=1,n; write(*,frmt) i,w(i),a(:,i); enddo
        deallocate(work)
        rn=0
        do i=0,n; do j=0,i; rn=rn+1; enddo; enddo
        rn=rn*2
        allocate(r(n,n),rho(rn),re(n,n))
!  convert density matrix
30      read(61,*,end=40) time, rho(:), psink
        ind=1
        r=0
        re=0
        pground=rho(1)
        do i=0,n; do j=0,i
         if(i.gt.0.and.j.gt.0)r(i,j)=rho(ind)+im*rho(ind+1)
         ind=ind+2;
        enddo; enddo
        do j=1,n; do i=1,j-1; r(i,j)=conjg(r(j,i)); enddo; enddo
        re=matmul(transpose(a),matmul(r,a))
        normalize=psink+pground
        do i=1,n; normalize=normalize+real(re(i,i)); enddo
        do i=1,n; do j=1,n; re(i,j)=re(i,j)/normalize; enddo; enddo
        write(*,*) time,pground,(real(re(i,i)),i=1,n),psink
        goto 30
40      write(*,*)
        end program eigenbasis
!gfortran -o eigenbasis eigenbasis.f -L /media/tom/My\ Passport/lapack/lapack-3.4.2 -l lapack -l refblas -O3 -fexternal-blas -ffixed-line-length-none
