!    A suite of routines to interpolate Model for Prediction Across
!    Scales (MPAS) Atmosphere (MPAS-A) core variables to a Gaussian
!    grid (of user specified spectral truncation) and readable by the
!    Gridded Analysis and Display System (GrADS).

!    Copyright (C) 2013 Henry R. Winterbottom

!    Email: Henry.Winterbottom@noaa.gov

!    Snail-mail:

!    Henry R. Winterbottom
!    NOAA/OAR/PSD R/PSD1
!    325 Broadway
!    Boulder, CO 80303-3328

!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation; either version 2 of
!    the License, or (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!    General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!    02110-1301 USA.

module mpi_interface

  !===================================================================

  !$$$ Module documentation block
  !

  !$$$

  !===================================================================

  use kinds

  !-------------------------------------------------------------------

  implicit none

  !-------------------------------------------------------------------

  ! Define necessary include files

  include "mpif.h"

  !-------------------------------------------------------------------

  ! Define all global MPI variables

  character                                              :: mpi_nodename(mpi_max_processor_name)
  character                                              :: mpi_noderequest
  logical                                                :: mpi_abort
  integer(kind=4),           dimension(:),   allocatable :: mpi_ranks
  integer(kind=4)                                        :: mpi_errorstatus(mpi_status_size)
  integer(kind=4)                                        :: mpi_masternode
  integer(kind=4)                                        :: mpi_slavenode
  integer(kind=4)                                        :: mpi_ierror
  integer(kind=4)                                        :: mpi_ierrorcode
  integer(kind=4)                                        :: mpi_procid
  integer(kind=4)                                        :: mpi_nprocs
  integer(kind=4)                                        :: mpi_node_source
  integer(kind=4)                                        :: mpi_node_destination
  integer(kind=4)                                        :: mpi_loopcount
  integer(kind=4)                                        :: mpi_request
  integer(kind=4)                                        :: mpi_group_user
  integer(kind=4)                                        :: mpi_group_nprocs
  integer(kind=4)                                        :: mpi_group_procid
  integer(kind=4)                                        :: mpi_group_begin
  integer(kind=4)                                        :: mpi_group_end

  !-------------------------------------------------------------------

contains

  !===================================================================
  
  ! mpi_interface_initialize.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_initialize()

    ! Initialize MPI session

    call mpi_init(mpi_ierror)

    ! Define rank for all nodes requested by user

    call mpi_comm_rank(mpi_comm_world,mpi_procid,mpi_ierror)

    ! Define the total number of nodes requested by user

    call mpi_comm_size(mpi_comm_world,mpi_nprocs,mpi_ierror)

    ! Define global variables

    mpi_masternode = 0

    ! Initialize global variables

    mpi_abort = .false.

  end subroutine mpi_interface_initialize

  !===================================================================

  ! mpi_interface_terminate.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_terminate()

    ! Terminate MPI session

    call mpi_finalize(mpi_ierror)
    
  end subroutine mpi_interface_terminate

  !===================================================================

  ! mpi_interface_define_comm.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_define_comm()

    ! Define variables computed within routine

    integer(kind=4),     dimension(:),   allocatable :: mpi_processes
    integer(kind=4)                                  :: mpi_worldgroup
    integer(kind=4)                                  :: mpi_newgroup

    ! Define counting variables

    integer                                          :: i, j, k
    integer                                          :: count

    !=================================================================

    ! Compute local variable

    mpi_group_nprocs = (mpi_group_end - mpi_group_begin) + 1

    ! Allocate memory for local variable

    if(.not. allocated(mpi_processes))                                 &
         & allocate(mpi_processes(mpi_group_nprocs))

    !-----------------------------------------------------------------

    ! Define local variable

    mpi_worldgroup = mpi_comm_world

    ! Initialize local variable

    count = 0

    ! Loop through each processor and define local variable

    do k = 1, mpi_nprocs

       ! Define local variable accordingly

       if(k .le. mpi_group_end .and. k .ge. mpi_group_begin) then

          ! Update local variable

          count = count + 1

          ! Define local variable

          mpi_processes(count) = k - 1

       end if ! if(k .le. mpi_group_end .and. k .ge. mpi_group_begin)

    end do ! do k = 1, mpi_nprocs

    ! Define local variables

    call mpi_comm_group(mpi_comm_world,mpi_worldgroup,mpi_ierror)
    call mpi_group_incl(mpi_worldgroup,mpi_group_nprocs,               &
         & mpi_processes(1:mpi_group_nprocs),mpi_newgroup,mpi_ierror)
    call mpi_comm_create(mpi_comm_world,mpi_newgroup,mpi_group_user,   &
         & mpi_ierror)

    !-----------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_processes)) deallocate(mpi_processes)

    !=================================================================

    ! Return calculated values

    return

    !=================================================================

  end subroutine mpi_interface_define_comm

  !===================================================================

end module mpi_interface
