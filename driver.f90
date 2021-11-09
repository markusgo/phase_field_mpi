program driver
	use parallel_pf
	implicit none
	integer :: ierror,nproc,procid,it,itmax
	real(8) :: dt,parsumphi,sumphi
	call mpi_start(ierror)
	call mpi_get_proc_info(ierror)

	call init_params(256,256,256,50.0d0,ierror)
	call init_phi(ierror)
	itmax=1000
	dt=1.0d-3
	parsumphi = 0.0d0
	sumphi = 0.0d0
	do it=1,itmax
		call euler_step_no_mech(dt)
		if(mod(it,50)==0) then
			parsumphi = sum(phi(:,:,1:zstride))
			call MPI_REDUCE(parsumphi,sumphi,1,MPI_REAL8,MPI_SUM,master,MPI_COMM_WORLD,ierror)
			if(ismaster) then
				print*,it,sumphi
			endif
		endif
	enddo
	call mpi_gather_grid(phi,phi_global)

	if(ismaster) then
		call save_vtk_scalar(phi_global,"phi","phi")
	endif
	call mpi_end(ierror)
end program