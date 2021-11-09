module parallel_pf
	implicit none
	include "mpif.h"

	integer :: lx,ly,lz
	real(8),allocatable,dimension(:,:,:) :: phi,dfdphi,phi_global
	real(8) :: radsph
	integer :: prank,psize,zstride
	integer,parameter :: master = 0

contains
	subroutine mpi_start(error)
		integer,intent(out) ::  error
		call MPI_INIT(error)
	end subroutine

	subroutine mpi_get_proc_info(error)
		integer,intent(out) :: error
		call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,error)
		call MPI_COMM_RANK(MPI_COMM_WORLD,prank,error)
		call MPI_BARRIER(MPI_COMM_WORLD,error)
		print*,"proc ",prank," of ",psize
	end subroutine

	subroutine mpi_end(error)
		integer,intent(out) ::  error
		call MPI_FINALIZE(error)
	end subroutine

	subroutine init_params(xdim,ydim,zdim,rsphere,error)
		integer,intent(in) :: xdim,ydim,zdim
		real(8),intent(in) :: rsphere
 		integer,intent(out) :: error
		lx = xdim
		ly = ydim
		lz = zdim
		radsph = rsphere
		call MPI_BARRIER(MPI_COMM_WORLD,error)
	end subroutine

	subroutine init_phi(error)
		integer,intent(out) :: error
		integer :: chunki,chunkf
		integer :: strides(psize),rcount(psize),i,j,k,z
		real(8) :: dx,dy,dz


		zstride = lz/psize
		chunki = prank*zstride+1
		chunkf = chunki+zstride
		

		allocate(phi(lx,ly,0:zstride+1))
		allocate(dfdphi(lx,ly,0:zstride+1))
		!do k=1,chunksize
		!	do j=1,ly
		!		do i=1,lx
		!			z = chunki+k
		!			dx = i-lx/2.0
		!			dy = j-ly/2.0
		!			dz = z-lz/2.0
		!			phi(i,j,k) = merge(1,-1,dx**2+dy**2+dz**2<radsph**2)
		!		enddo
		!	enddo
		!enddo
		phi = prank
		
		if(prank==master) then
			allocate(phi_global(lx,ly,lz))
			do i=1,psize
				strides(i) = (i-1)*lx*ly*zstride
				rcount(i) = lx*ly*zstride
			enddo
		endif
		call MPI_BARRIER(MPI_COMM_WORLD,error)
		
		call MPI_GATHERV(phi(:,:,1:zstride),lx*ly*zstride,MPI_REAL8, &
						 phi_global,rcount,strides,MPI_REAL8, &
						 master,MPI_COMM_WORLD,error)
	end subroutine

	subroutine ghostswap(data_swap)
		real(8),intent(inout) :: data_swap(lx,ly,0:zstride+1) 
		integer :: vizprev,viznext,error,datasize
		datasize = lx*ly
		viznext = merge(0,prank+1,prank==psize-1)
		vizprev = merge(psize-1,prank-1,prank==0)

		call MPI_SENDRECV(data_swap(:,:,1)   ,datasize,MPI_REAL8,vizprev,1, &
						  data_swap(:,:,zstride+1),datasize,MPI_REAL8,viznext,1, &
						  MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
		call MPI_SENDRECV(data_swap(:,:,zstride),datasize,MPI_REAL8,viznext,1, &
						  data_swap(:,:,0) ,datasize,MPI_REAL8,vizprev,1, &
						  MPI_COMM_WORLD,MPI_STATUS_IGNORE,error)
	end subroutine

	subroutine euler_step_no_mech(dt)
		real(dp),intent(in) :: dt
		integer :: i,j,k,z
		real(dp) :: lapphi,lapdfdphi

		do k = 1,lz
			do j = 1,ly
				do i = 1,lx
					z = 
					lapphi = phi(inext,j,k) + phi(iprev,j,k) &
						   + phi(i,jnext,k) + phi(i,jprev,k) &
						   + phi(i,j,knext) + phi(i,j,kprev) &
						   - sixdp * phi(i,j,k)
					dfdphi(i,j,k) = rho_phi * ( -phi(i,j,k) + phi(i,j,k)**3 - epsilon**2*lapphi)
				enddo
			enddo
		enddo

		do k = 1,lz
			do j = 1,ly
				do i = 1,lx
					inext = modulo(i,lx) + 1
					iprev = modulo(i-2,lx) + 1
					jnext = modulo(j,ly) + 1
					jprev = modulo(j-2,ly) + 1
					knext = modulo(k,lz) + 1
					kprev = modulo(k-2,lz) + 1
					lapdfdphi = dfdphi(inext,j,k) + dfdphi(iprev,j,k) &
						      + dfdphi(i,jnext,k) + dfdphi(i,jprev,k) &
						      + dfdphi(i,j,knext) + dfdphi(i,j,kprev) &
						      - sixdp * dfdphi(i,j,k)
					phinew(i,j,k) = phi(i,j,k) + dt * mob * lapdfdphi
				enddo
			enddo
		enddo
		phi = phinew
	end subroutine

	subroutine save_vtk_scalar(mat,name,dataset)
    	implicit none
		integer :: i,j,k
		real(8),intent(in) :: mat(lx,ly,lz)
		character(len=*),intent(in) :: name,dataset

		open(unit=13,file=trim(name)//'.vtk')
		write(13,"(A)") "# vtk DataFile Version 2.0"
		write(13,"(A)") "test file"
		write(13,"(A)") "ASCII"
		write(13,"(A)") "DATASET STRUCTURED_POINTS"
		write(13,"(A,1x,I0,1x,I0,1x,I0)") "DIMENSIONS",lx,ly,lz
		write(13,"(A)") "ORIGIN 0.0 0.0 0.0"
		write(13,"(A)") "SPACING 1.0 1.0 1.0"
		write(13,"(A)") ""
		write(13,"(A,1x,I0)") "POINT_DATA",lx*ly*lz
		write(13,"(A,A,A)") "SCALARS ",trim(dataset)," float"
		write(13,"(A)") "LOOKUP_TABLE default"
	    do k=1,lz
	    	do j=1,ly
	    		do i=1,lx
	    			write(13,'(F8.4,1X)',advance='no') mat(i,j,k)
	    		enddo
	    	enddo
	    enddo
    	close(13)
	end subroutine
end module