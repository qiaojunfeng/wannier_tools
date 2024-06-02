subroutine readNormalSpinR()
   !>> Read in the tight-binding model for spin operator from wannier90.
   ! This file can be obtained by transforming the gauge of w90 spn file
   ! from Bloch gauge to Wannier gauge, and write the x,y,z components
   ! to three files, in the wannier90 _hr.dat format.
   !
   ! License: GPL V3

   use para
   !> in: N of wann
   !> out : nth atom

   implicit none

   character*4 :: c_temp

   integer :: i, j, ir, ia, io 
   integer :: i1, i2, i3, i4, i5
   integer :: n, m, ir0
   integer :: add_electric_field
   integer :: nwann, nwann_nsoc

   real(dp) :: static_potential
   real(dp) :: tot, rh, ih
   real(dp) :: pos(Origin_cell%Num_atoms)

   character(80) :: Spinfile  = "spin_hr.dat"           ! filename
   character(1) :: ijk2xyz(3) = ['x', 'y', 'z']
   integer, allocatable     :: tmp_ndegen(:), tmp_irvec(:,:)
   integer :: ispin

   !> add a check for NumOccupied parameter
   if (NumOccupied<=0 .or. NumOccupied>Num_wann) then
      write(stdout, '(a, i6, a)')">>> ERROR: NumOccupied should be in [1, ", Num_wann, " ]"
      write(stdout, '(a)')">>> Usually, it is the number of occupied Wannier bands."
      stop
   endif

    do ispin = 1, 3
        if(cpuid.eq.0) write(stdout,*)' Read spin operator ', ijk2xyz(ispin)
        open(12, file=trim(Spinfile) // "." // ijk2xyz(ispin), status='OLD')
        !> for Normal HmnR obtained from Wannier90 or sparse HmnR

        !> skip a comment line
        read(12, *)

        !> number of Wannier orbitals in the hr file
        nwann=0
        read(12, *)nwann
        if (nwann==0) then
            stop "ERROR : num_wann is zero in hr file"
        endif
        nwann_nsoc=nwann
        if (SOC>0) nwann_nsoc= nwann/2

        if ((soc==0 .and. sum(Origin_cell%nprojs)/=nwann .and. .not.Add_Zeeman_Field) .or. &
        (soc>0 .and. sum(Origin_cell%nprojs)/=nwann/2))then
        print *, 'sum(Origin_cell%nprojs), num_wann, num_wann/2'
        print *, sum(Origin_cell%nprojs), nwann, nwann/2
        print *, "ERROR: Maybe the SOC tags in the SYSTEM is wrongly set"
        stop "ERROR: the summation of all projectors times spin degeneracy is not equal to num_wann"
        endif

        !> number of lattice vectors taken into account
        read(12, *) Nrpts
        if (.not. allocated(SpinR)) allocate(SpinR(Num_wann,Num_wann,nrpts,3))
        if (.not. allocated(tmp_irvec)) allocate(tmp_irvec(3,nrpts))
        if (.not. allocated(tmp_ndegen)) allocate(tmp_ndegen(nrpts))
        tmp_irvec = 0
        tmp_ndegen = 1

        !> The degeneracy of each R point, this is caused by the Fourier transform in the Wigner-Seitz cell
        read(12,*,end=1001) (tmp_ndegen(i), i=1, Nrpts)
        do ir = 1, nrpts
            if (tmp_ndegen(ir) /= ndegen(ir)) then
                stop "ERROR : ndegen is different from Hamiltonian hr.dat file"
            endif
        enddo
        do ir=1,Nrpts
            do n=1,nwann
                do m=1,nwann
                    read(12,*,end=1001) i1, i2, i3, i4, i5, rh, ih
                    if ((irvec(1,ir) /= i1) .or. &
                        (irvec(2,ir) /= i2) .or. &
                        (irvec(3,ir) /= i3) .or. &
                        (i4 /= m) .or. (i5 /= n)) then
                        stop "ERROR : irvec is different from Hamiltonian hr.dat file"
                    ENDIF
                    SpinR(i4,i5,ir,ispin) = dcmplx(rh,ih)
                end do
            enddo
        enddo

        1001 continue
        close(12)
    enddo

   !check sum rule
   tot= 0d0
   do ir=1, Nrpts
      tot= tot+ 1d0/ndegen(ir)
   enddo

   call get_SpinR_cell(Cell_defined_by_surface)

   return
end subroutine readNormalSpinR

subroutine get_SpinR_cell(cell)
    !> Get new SpinR for a new cell with the same size as the previous one
    ! must be called after get_hmnr_cell
    use para
    implicit none
 
    type(cell_type) :: cell
 
    !type(dense_tb_hr) :: cell_hr
 
    integer :: ir, i, j, iter
    real(dp) :: shift_vec_direct(3)
 
    !> for newcell
    real(dp) :: apos1d(3),apos2d(3),apos1dprime(3),apos2dprime(3)
    !>count newcell nrpts
    integer :: max_ir
    integer :: nir1,nir2,nir3, ir_cell
    integer :: nrpts_new, nrpts_max
    real(dp) :: new_ia, new_ib, new_ic, max_val
 
    !> all new irs
    integer, allocatable  :: rpts_array(:, :, :), rpts_map(:, :, :)
    integer :: irn1(3),irn2(3)
 
    !> The Allocate Process
    integer :: ia1,ia2,ia1prime,ia2prime
 
    max_ir=8
    nrpts_max=(2*max_ir+1)**3
    allocate( rpts_array(-max_ir:max_ir,-max_ir:max_ir,-max_ir:max_ir))
    allocate( rpts_map(-max_ir:max_ir,-max_ir:max_ir,-max_ir:max_ir))
 
    ! allocate(Rmn_old(3))
    ! allocate(Rmn_new(3))
    ! allocate(irvec_new(3))
    ! allocate(irvec_new_int(3))
 
    ! call cart_direct_real(shift_to_topsurface_cart, shift_vec_direct, cell%lattice)
 
    call date_and_time(DATE=date_now,ZONE=zone_now, TIME=time_now)
    !> Get new Hrs
    rpts_array=0
    ! nrpts_surfacecell=0
    rpts_map= 0
 
    !> 1. Get number of R points for the new cell
    do ir=1, nrpts
       do j=1, Num_wann
          do i=1, Num_wann
             ia1 = Origin_cell%spinorbital_to_atom_index(i)
             ia2 = Origin_cell%spinorbital_to_atom_index(j)
             apos1d = Origin_cell%Atom_position_direct(:, ia1)
             apos2d = Origin_cell%Atom_position_direct(:, ia2)
 
             ia1prime = Cell_defined_by_surface%spinorbital_to_atom_index(i)
             ia2prime = Cell_defined_by_surface%spinorbital_to_atom_index(j)
             apos1dprime = Cell_defined_by_surface%Atom_position_direct(:,ia1prime)
             apos2dprime = Cell_defined_by_surface%Atom_position_direct(:,ia2prime)
 
             Rmn_old = irvec(:, ir) + apos2d - apos1d
             call latticetransform(Rmn_old(1),Rmn_old(2),Rmn_old(3),Rmn_new(1),Rmn_new(2),Rmn_new(3))
             irvec_new = Rmn_new - (apos2dprime - apos1dprime)
 
             !> Due to the accuracy of computing, need to rounding  (but always tiny)
             irvec_new_int(1) = ANINT(irvec_new(1))
             irvec_new_int(2) = ANINT(irvec_new(2))
             irvec_new_int(3) = ANINT(irvec_new(3))
 
             if (abs(irvec_new_int(1))>max_ir .or. abs(irvec_new_int(2))>max_ir .or. abs(irvec_new_int(3))>max_ir) cycle
             rpts_array(irvec_new_int(1),irvec_new_int(2),irvec_new_int(3))=1
 
          enddo
       enddo
    enddo
 
    !> The total number of lattice points searched above
    if (nrpts_surfacecell /= sum(rpts_array)) then
        stop "ERROR : nrpts_surfacecell is different"
    endif
 
    !> 2. Create an order map to sign the R points generated in step1
    iter = 0
    do nir3=-max_ir,max_ir
       do nir2=-max_ir,max_ir
          do nir1=-max_ir,max_ir
             if (rpts_array(nir1, nir2, nir3)==1) then
                iter=iter+1
                rpts_map(nir1, nir2, nir3)=iter
             endif
          enddo
       enddo
    enddo

    allocate(SpinR_surfacecell(Num_wann, Num_wann, nrpts_surfacecell, 3))
    ! allocate(ndegen_surfacecell(nrpts_surfacecell))
    SpinR_surfacecell= 0d0
    ! ndegen_surfacecell= 1
 
    !> 3. Allocate the old HmnR to the new HmnR.
    !>  Note: The cell can't be treated as a mass point. We need to consider the relative coordinates of atoms.
    do ir=1, nrpts
       do j=1, Num_wann
          do i=1, Num_wann
 
             ia1 = Origin_cell%spinorbital_to_atom_index(i)
             ia2 = Origin_cell%spinorbital_to_atom_index(j)
             apos1d = Origin_cell%Atom_position_direct(:, ia1)
             apos2d = Origin_cell%Atom_position_direct(:, ia2) 
 
             ia1prime = Cell_defined_by_surface%spinorbital_to_atom_index(i)
             ia2prime = Cell_defined_by_surface%spinorbital_to_atom_index(j)
             apos1dprime = Cell_defined_by_surface%Atom_position_direct(:,ia1prime)
             apos2dprime = Cell_defined_by_surface%Atom_position_direct(:,ia2prime)
 
             !> R'_mn = R' + tau'_2 - tau'_1
             !> R_mn = R + tau_2 - tau_1
             !> R'_mn = Pinv * R_mn
             !> R' = (Pinv * R_mn) - (tau'_2 - tau'_1)
 
             Rmn_old = irvec(:, ir) + apos2d - apos1d
             call latticetransform(Rmn_old(1),Rmn_old(2),Rmn_old(3),Rmn_new(1),Rmn_new(2),Rmn_new(3))
             irvec_new = Rmn_new - (apos2dprime - apos1dprime)
 
             !> For safety, we perform a rounding operation due to the accuracy of computing
             irvec_new_int(1) = ANINT(irvec_new(1))
             irvec_new_int(2) = ANINT(irvec_new(2))
             irvec_new_int(3) = ANINT(irvec_new(3))
 
             if (abs(irvec_new_int(1))>max_ir .or. abs(irvec_new_int(2))>max_ir .or. abs(irvec_new_int(3))>max_ir) cycle
             ir_cell = rpts_map(irvec_new_int(1), irvec_new_int(2), irvec_new_int(3))
             SpinR_surfacecell(i,j,ir_cell,:) = SpinR(i,j,ir,:)/ndegen(ir)
 
          enddo
       enddo
    enddo
 
    return
 end subroutine get_SpinR_cell
