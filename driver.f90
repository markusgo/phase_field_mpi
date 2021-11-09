program driver
	use parallel_pf
	implicit none
	integer :: ierror,nproc,procid

	call mpi_start(ierror)
	call mpi_get_proc_info(ierror)

	call init_params(4,4,4,20.0d0,ierror)
	call init_phi(ierror)
	call ghostswap(phi)

	!if(prank==master) then
	!	call save_vtk_scalar(phi_global,"phi","phi")
	!endif

	if(prank==0) then
		print*,phi(:,:,2)
	endif

	call mpi_end(ierror)
end program