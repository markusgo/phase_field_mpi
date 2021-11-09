program driver
	use parallel_pf
	implicit none
	integer :: ierror,nproc,procid,it,itmax
	real(8) :: dt

	call mpi_start(ierror)
	call mpi_get_proc_info(ierror)

	call init_params(256,256,256,50.0d0,ierror)
	call init_phi(ierror)
	itmax=1000
	dt=1.0d-3
	do it=1,itmax
		call euler_step_no_mech(dt)
	enddo
	call mpi_gather_grid(phi,phi_global)

	if(prank==master) then
		call save_vtk_scalar(phi_global,"phi","phi")
	endif
	call mpi_end(ierror)
end program