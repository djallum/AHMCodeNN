module Test
  implicit none
  integer :: dim = 10

contains
  subroutine DefineCluster( SitePotential, Sites, weakL, weakR, ClusterSize )
    implicit none

    ! Inputs
    real, dimension(dim), intent(in) :: SitePotential
    integer, intent(in) :: weakL, weakR
    integer, intent(in) :: ClusterSize

    ! Outputs
    real, dimension(ClusterSize), intent(out) :: Sites

    ! Other
    integer :: EndSite, FirstSite
       
    
    if ( weakL .gt. weakR ) then
       EndSite = dim - weakL
       Sites(1:EndSite) = SitePotential((weakL):(dim - sitesremoved))
       Sites((EndSite+1):ClusterSize) = SitePotential(1:weakR)
    else if ( (weakL .gt. weakR) .and. ( ClusterSize .eq. 1 ) ) then
       Sites = SitePotential(weakR:weakR)
    else if ( size(SitePotential) .eq. ClusterSize ) then
       Sites = SitePotential
    else
       Sites = SitePotential((weakL+1):weakR)
    end if

  end subroutine DefineCluster

end module Test

program Prog
  implicit none

  integer :: i, weakL, weakR, ClusterSize
  real :: SitePotential(dim), Sites(4)

  weakL = 7
  
  
  do i=1,dim
     SitePotential(i) = i
  end do

  
