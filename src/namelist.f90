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

module namelist

  !=======================================================================

  use kinds

  !-----------------------------------------------------------------------

  use mpi_interface

  !-----------------------------------------------------------------------

  implicit none
  
  !-----------------------------------------------------------------------

  ! Define global variables

  character(len=500)                             :: variable_filename         = 'NOT USED'
  logical                                        :: debug                     = .false.
  logical                                        :: is_cori                   = .false.
  logical                                        :: is_psfc                   = .false.
  logical                                        :: is_pslp                   = .false.
  logical                                        :: is_u10m                   = .false.
  logical                                        :: is_v10m                   = .false.
  logical                                        :: is_t2m                    = .false.
  logical                                        :: is_sst                    = .false.
  logical                                        :: is_ter                    = .false.
  logical                                        :: is_uwnd                   = .false.
  logical                                        :: is_vwnd                   = .false.
  logical                                        :: is_vort                   = .false.
  logical                                        :: is_divg                   = .false.
  logical                                        :: is_rh                     = .false.
  logical                                        :: is_qv                     = .false.
  logical                                        :: is_temp                   = .false.
  logical                                        :: is_vtemp                  = .false.
  logical                                        :: is_theta                  = .false.
  logical                                        :: is_thetae                 = .false.
  logical                                        :: is_psi                    = .false.
  logical                                        :: is_chi                    = .false.
  logical                                        :: is_pwat                   = .false.
  logical                                        :: is_pv                     = .false.
  logical                                        :: is_dbz                    = .false.
  logical                                        :: is_cdbz                   = .false.
  logical                                        :: is_ght                    = .false.
  logical                                        :: is_pinterp                = .false.
  logical                                        :: is_zinterp                = .false.
  logical                                        :: is_teinterp               = .false.
  logical                                        :: is_smooth                 = .false.
  real(r_kind)                                   :: plevs(100)                = -999.0
  real(r_kind)                                   :: zlevs(100)                = -999.0
  real(r_kind)                                   :: televs(100)               = -999.0
  real(r_kind)                                   :: barnes_weights_threshold  = 1.e-5
  real(r_kind)                                   :: barnes_distance_threshold = 300000.0
  integer                                        :: barnes_nneighbors         = 10
  integer                                        :: barnes_npasses            = 2
  integer                                        :: grid_mtrunc               = 62
  integer                                        :: nlev                      = 0
  integer                                        :: namelist_io_error         = 0


  namelist /fileio/ debug, variable_filename
  namelist /variableio/ is_psfc, is_pslp, is_u10m, is_v10m, is_t2m,      &
       & is_sst, is_ter, is_uwnd, is_vwnd, is_rh, is_vtemp, is_temp,     &
       & is_theta, is_thetae, is_vort, is_divg, is_psi, is_chi, is_pwat, &
       & is_pv, is_qv, is_cori, is_dbz, is_cdbz, is_ght, is_smooth
  namelist /interpio/ grid_mtrunc, barnes_nneighbors,                    &
       & barnes_npasses, barnes_distance_threshold,                      &
       & barnes_weights_threshold
  namelist /pinterpio/ is_pinterp, plevs
  namelist /zinterpio/ is_zinterp, zlevs
  namelist /teinterpio/ is_teinterp, televs

  !---------------------------------------------------------------------

contains

  !=======================================================================

  ! namelistparams.f90:

  !-----------------------------------------------------------------------

  subroutine namelistparams()
    
    ! Define variables computed within routine

    logical                                        :: is_it_there
    integer                                        :: unit_nml

    ! Define counting variables

    integer                                        :: i, j, k
    
    !=====================================================================

    ! Initialize local variables

    unit_nml    = 9
    is_it_there = .false.

    ! Define local variable

    inquire(file='mpastograds.input',exist = is_it_there)

    ! Check local variable and proceed accordingly

    if(is_it_there) then

       ! Open external file

       open(file = 'mpastograds.input',                                   &
            unit = unit_nml        ,                                      &
            status = 'old'         ,                                      &
            form = 'formatted'     ,                                      &
            action = 'read'        ,                                      &
            access = 'sequential'  )

       ! Define local variables

       read(unit_nml,NML = fileio)
       read(unit_nml,NML = variableio)
       read(unit_nml,NML = interpio)
       read(unit_nml,NML = pinterpio)
       read(unit_nml,NML = zinterpio)
       read(unit_nml,NML = teinterpio)
       
       ! Close external file
    
       close(unit_nml)

       ! Compute local variable accordingly

       if(is_pinterp) then

          ! Initialize local variable

          nlev = 0

          ! Loop through all user specified vertical levels

          do k = 1, 100

             ! Update local variable accordingly

             if(plevs(k) .ne. -999.0) then

                ! Rescale local variables

                plevs(k) = plevs(k)*100.0
                
                ! Update local variable

                nlev = nlev + 1
          
             end if ! if(plevs(k) .ne. -999.0)

          end do ! do k = 1, 100

       end if ! if(is_pinterp)

       ! Compute local variable accordingly

       if(is_zinterp) then

          ! Initialize local variable

          nlev = 0

          ! Loop through all user specified vertical levels

          do k = 1, 100

             ! Update local variable accordingly

             if(zlevs(k) .ne. -999.0) then
                
                ! Update local variable

                nlev = nlev + 1
          
             end if ! if(zlevs(k) .ne. -999.0)

          end do ! do k = 1, 100

       end if ! if(is_zinterp)

       ! Compute local variable accordingly

       if(is_teinterp) then

          ! Initialize local variable

          nlev = 0

          ! Loop through all user specified vertical levels

          do k = 1, 100

             ! Update local variable accordingly

             if(televs(k) .ne. -999.0) then
                
                ! Update local variable

                nlev = nlev + 1
          
             end if ! if(televs(k) .ne. -999.0)

          end do ! do k = 1, 100

       end if ! if(is_teinterp)

    end if ! if(is_it_there)

    ! Check local variable and proceed accordingly

    if(.not. is_it_there) then 

       ! Print message to user

       write(6,500)

       ! Update local variable

       namelist_io_error = 1.0

    end if ! if(.not. is_it_there)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Print message to user

       write(6,*) '&FILEIO'
       write(6,*) 'DEBUG             = ', debug
       write(6,*) 'VARIABLE_FILENAME = ',                                 &
            & trim(adjustl(variable_filename))
       write(6,*) '/'
       write(6,*) '&VARIABLEIO'
       if(is_cori)   write(6,*) 'IS_CORI   = ', is_cori
       if(is_psfc)   write(6,*) 'IS_PSFC   = ', is_psfc
       if(is_pslp)   write(6,*) 'IS_PSLP   = ', is_pslp
       if(is_u10m)   write(6,*) 'IS_U10M   = ', is_u10m
       if(is_v10m)   write(6,*) 'IS_V10M   = ', is_v10m
       if(is_t2m)    write(6,*) 'IS_T2M    = ', is_t2m
       if(is_sst)    write(6,*) 'IS_SST    = ', is_sst
       if(is_ter)    write(6,*) 'IS_TER    = ', is_ter
       if(is_uwnd)   write(6,*) 'IS_UWND   = ', is_uwnd
       if(is_vwnd)   write(6,*) 'IS_VWND   = ', is_vwnd
       if(is_vort)   write(6,*) 'IS_VORT   = ', is_vort
       if(is_divg)   write(6,*) 'IS_DIVG   = ', is_divg
       if(is_rh)     write(6,*) 'IS_RH     = ', is_rh
       if(is_qv)     write(6,*) 'IS_QV     = ', is_qv
       if(is_temp)   write(6,*) 'IS_TEMP   = ', is_temp
       if(is_vtemp)  write(6,*) 'IS_VTEMP  = ', is_vtemp
       if(is_theta)  write(6,*) 'IS_THETA  = ', is_theta
       if(is_thetae) write(6,*) 'IS_THETAE = ', is_thetae
       if(is_psi)    write(6,*) 'IS_PSI    = ', is_psi
       if(is_chi)    write(6,*) 'IS_CHI    = ', is_chi
       if(is_pwat)   write(6,*) 'IS_PWAT   = ', is_pwat
       if(is_pv)     write(6,*) 'IS_PV     = ', is_pv
       if(is_dbz)    write(6,*) 'IS_DBZ    = ', is_dbz
       if(is_cdbz)   write(6,*) 'IS_CDBZ   = ', is_cdbz
       if(is_ght)    write(6,*) 'IS_GHT    = ', is_ght
       write(6,*) 'IS_SMOOTH = ', is_smooth
       write(6,*) '/'
       write(6,*) '&INTERPIO'
       write(6,*) 'GRID_MTRUNC               = ', grid_mtrunc
       write(6,*) 'BARNES_NNEIGHBORS         = ', barnes_nneighbors
       write(6,*) 'BARNES_NPASSES            = ', barnes_npasses
       write(6,*) 'BARNES_DISTANCE_THRESHOLD = ',                         &
            & barnes_distance_threshold
       write(6,*) 'BARNES_WEIGHTS_THRESHOLD  = ',                         &
            & barnes_weights_threshold
       write(6,*) '/'
       if(is_pinterp) then
          write(6,*) '&PINTERPIO'
          write(6,*) 'IS_PINTERP = ', is_pinterp
          write(6,*) 'PLEVS      = ', plevs(1:nlev)
          write(6,*) '/'
       end if ! if(is_pinterp)
       if(is_zinterp) then
          write(6,*) '&ZINTERPIO'
          write(6,*) 'IS_ZINTERP = ', is_zinterp
          write(6,*) 'ZLEVS      = ', zlevs(1:nlev)
          write(6,*) '/'
       end if ! if(is_zinterp)
       if(is_teinterp) then
          write(6,*) '&TEINTERPIO'
          write(6,*) 'IS_TEINTERP = ', is_teinterp
          write(6,*) 'TELEVS      = ', televs(1:nlev)
          write(6,*) '/'
       end if ! if(is_teinterp)
       write(6,*) ' '

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations
    
    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Define format statements

500 format('NAMELISTPARAMS: mpas-post.input not found in current ',        &
         & 'working directory. ABORTING!!!!')

    !=====================================================================

  end subroutine namelistparams
  
  !=======================================================================

end module namelist
