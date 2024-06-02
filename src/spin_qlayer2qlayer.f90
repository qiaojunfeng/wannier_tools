  subroutine spin_qlayer2qlayer(k, sigma_x, sigma_y, sigma_z)
     ! This subroutine caculates spin operator between
     ! slabs For surface state calculation

     use para

     implicit none

     ! loop index
     integer :: i,j,iR

     ! index used to sign irvec     
     real(dp) :: ia,ib,ic

     ! new index used to sign irvec     
     real(dp) :: new_ia,new_ib,new_ic
     integer :: inew_ic
     integer :: int_ic

     ! wave vector k times lattice vector R  
     real(Dp) :: kdotr  

     ! input wave vector k's cooridinates
     real(Dp),intent(in) :: k(2)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin

     complex(dp) :: ratio

     complex(Dp), allocatable :: Sij(:, :, :, :)

     ! H00 Hamiltonian between nearest neighbour-quintuple-layers
     ! the factor 2 is induced by spin

     !     complex(Dp),allocatable,intent(out) :: H00new(:,:)
     complex(Dp),intent(out) :: sigma_x(Ndim,Ndim)
     complex(Dp),intent(out) :: sigma_y(Ndim,Ndim)
     complex(Dp),intent(out) :: sigma_z(Ndim,Ndim)

     allocate(Sij(Num_wann,Num_wann,-ijmax:ijmax, 3))

     Sij=0.0d0
     do iR=1,nrpts_surfacecell
        ia=irvec_surfacecell(1,iR)
        ib=irvec_surfacecell(2,iR)
        ic=irvec_surfacecell(3,iR)

        int_ic= int(ic)
        if (abs(ic).le.ijmax)then
           kdotr=k(1)*ia+ k(2)*ib
           ratio=cos(2d0*pi*kdotr)+zi*sin(2d0*pi*kdotr)

           Sij(1:Num_wann, 1:Num_wann, int_ic, :)&
           =Sij(1:Num_wann, 1:Num_wann, int_ic, :)&
           +SpinR_surfacecell(:,:,iR,:)*ratio/ndegen_surfacecell(iR)
        endif

     enddo

     sigma_x=0.0d0
     sigma_y=0.0d0
     sigma_z=0.0d0

     ! nslab's principle layer 
     do i=1,Np
     do j=1,Np
        if (abs(i-j).le.(ijmax)) then
          sigma_x(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Sij(:,:,j-i,1)
          sigma_y(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Sij(:,:,j-i,2)
          sigma_z(Num_wann*(i-1)+1:Num_wann*i,Num_wann*(j-1)+1:Num_wann*j)&
                =Sij(:,:,j-i,3)
        endif
     enddo
     enddo

     deallocate(Sij)

  return
  end subroutine spin_qlayer2qlayer
