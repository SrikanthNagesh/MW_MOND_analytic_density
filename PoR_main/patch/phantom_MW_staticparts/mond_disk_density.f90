!########################################################
!########################################################
!########################################################
!########################################################

!!!! This file is added by STN !!!!!!!

!========================================================
module mond_disk_density

  use amr_parameters
  use mond_disk_parameters
  use mond_parameters

  implicit none

  real(dp) :: G_kpc
  real(dp) :: a0_code

  real(dp) :: M_inner
  real(dp) :: M_outer

  integer, parameter :: nr_table=20000

  real(dp), dimension(nr_table) :: r_table
  real(dp), dimension(nr_table) :: logrho_table

  real(dp) :: log_rmin
  real(dp) :: dlogr

  logical :: disk_table_ready=.false.

contains

!========================================================
subroutine init_mond_disk(mass_disk,f_inner,f_outer)

  implicit none

  real(dp),intent(in)::mass_disk
  real(dp),intent(in)::f_inner,f_outer

  real(dp)::scale_l,scale_t
  real(dp)::scale_d,scale_v
  real(dp)::scale_nH,scale_T2

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! RAMSES code units: G = 1
  G_kpc = 1.0d0

  ! MOND acceleration scale in code units
  if (a0_ms2 > 0.d0) then
     a0_code = a0_ms2*100.d0/scale_l*scale_t**2
  else
     a0_code = a0_mond*1.d3
  endif

  M_inner = f_inner*mass_disk
  M_outer = f_outer*mass_disk

  write(*,*) "MOND initialization"
  write(*,*) "G_kpc   =",G_kpc
  write(*,*) "a0_code =",a0_code
  write(*,*) "M_inner =",M_inner/1.0d10,"10^10 Msun"
  write(*,*) "M_outer =",M_outer/1.0d10,"10^10 Msun"

end subroutine init_mond_disk


!========================================================
! Construct TABLE OF THE PHYSICAL BARYONIC DENSITY
!
! NO MOND correction is applied here.
!
! rho(r) =
!
! M_inner /(4*pi*Rd_inner^2*r) * exp(-r/Rd_inner)
!
! +
!
! M_outer /(4*pi*Rd_outer^2*r) * exp(-r/Rd_outer)
!
!========================================================
subroutine build_density_table()

  implicit none

  integer :: i
  real(dp) :: r
  real(dp) :: rho1,rho2,rho
  real(dp), parameter :: rmin=1.d-5

  log_rmin = log(rmin)
  dlogr = (log(rmax)-log(rmin))/real(nr_table-1,dp)

  do i=1,nr_table

     r = exp(log_rmin + real(i-1,dp)*dlogr)

     r_table(i)=r

     call rho_disk_component(M_inner,r,Rd_inner,rho1)
     call rho_disk_component(M_outer,r,Rd_outer,rho2)

     rho = rho1 + rho2

     if (rho > 0.d0) then
        logrho_table(i)=log(rho)
     else
        logrho_table(i)=log(1.d-50)
     endif

  enddo

  disk_table_ready=.true.

  write(*,*) "BARYONIC density table constructed"
  write(*,*) "r range:",r_table(1),r_table(nr_table)
  write(*,*) "rho range:",exp(logrho_table(1)), &
                         exp(logrho_table(nr_table))

  open(100,file='mond_density_table.dat',status='replace')

  do i=1,nr_table
     write(100,*) r_table(i), exp(logrho_table(i))
  enddo

  close(100)

end subroutine build_density_table


!========================================================
! Enclosed mass fraction
!
! M(<r) = M * f(r,rd)
!
!========================================================
subroutine f_disk(r,rd,f)

  implicit none

  real(dp),intent(in)::r,rd
  real(dp),intent(out)::f

  f = 1.d0-exp(-r/rd)*(1.d0+r/rd)

end subroutine f_disk


!========================================================
! Derivative of enclosed mass fraction
!
! df/dr = r/rd^2 * exp(-r/rd)
!
!========================================================
subroutine df_dr_disk(r,rd,df)

  implicit none

  real(dp),intent(in)::r,rd
  real(dp),intent(out)::df

  df = r/(rd*rd)*exp(-r/rd)

end subroutine df_dr_disk


!========================================================
! PHYSICAL BARYONIC DENSITY OF ONE COMPONENT
!
! dM/dr = M * df/dr
!
! rho = 1/(4*pi*r^2) * dM/dr
!
! This simplifies to
!
! rho = M/(4*pi*rd^2*r) * exp(-r/rd)
!
! NO MOND CORRECTION.
!
!========================================================
subroutine rho_disk_component(M,r,rd,rho)

  implicit none

  real(dp),intent(in)::M,r,rd
  real(dp),intent(out)::rho

  real(dp)::pi

  pi=acos(-1.d0)

  rho = M/(4.d0*pi*rd*rd*r)*exp(-r/rd)

end subroutine rho_disk_component


!========================================================
! TOTAL PHYSICAL BARYONIC DENSITY
!
! This is the quantity to give to the QUMOND solver.
!
!========================================================
subroutine rho_disk_ana(r,rho)

  implicit none

  real(dp),intent(in)::r
  real(dp),intent(out)::rho

  real(dp)::rho1,rho2

  call rho_disk_component(M_inner,r,Rd_inner,rho1)

  call rho_disk_component(M_outer,r,Rd_outer,rho2)

  rho = rho1 + rho2

end subroutine rho_disk_ana


!========================================================
! Interpolate baryonic density table
!========================================================
subroutine interpolate_density(r,rho)

  implicit none

  real(dp),intent(in) :: r
  real(dp),intent(out):: rho

  integer :: i
  real(dp):: x
  real(dp):: y0,y1,y2,y3
  real(dp):: logr
  real(dp):: logrho_table_value

  if(.not.disk_table_ready) then
     write(*,*) 'ERROR: MOND density table not initialized'
     call clean_stop
  endif

  if(r >= rmax) then
     rho=1.d-5
     return
  endif

  logr=log(r)

  x=(logr-log_rmin)/dlogr

  i=int(x)+1

  if(i<2) i=2
  if(i>nr_table-2) i=nr_table-2

  x=x-real(i-1,dp)

  y0=logrho_table(i-1)
  y1=logrho_table(i)
  y2=logrho_table(i+1)
  y3=logrho_table(i+2)

  logrho_table_value = y1 + 0.5d0*x*( &
       y2-y0 + x*( &
       2.d0*y0-5.d0*y1+4.d0*y2-y3 + x*( &
       3.d0*(y1-y2)+y3-y0 )))

  rho=exp(logrho_table_value)

end subroutine interpolate_density


!========================================================
! Print ANALYTICAL MOND rotation curve for comparison.
!
! This does NOT modify the density.
!
!========================================================
subroutine print_mond_rotation_table()

  use amr_parameters
  use mond_disk_parameters

  implicit none

  integer :: i,nr
  real(dp) :: r
  real(dp) :: f1,f2
  real(dp) :: mb_inner,mb_outer,mb_total
  real(dp) :: gN,y,nu
  real(dp) :: mc
  real(dp) :: vc,vc_newt
  real(dp) :: rmin

  rmin = 1.0d-2

  nr = 200

  open(unit=99,file='mond_rotation_table.dat',status='replace')

  write(99,'(A)') &
  '# r(kpc) Mb_inner Mb_outer Mb_total gN y nu Mc vc vc_newt'

  do i=1,nr

     r = rmin * (rmax/rmin)**((i-1.d0)/(nr-1.d0))

     call f_disk(r,Rd_inner,f1)
     call f_disk(r,Rd_outer,f2)

     mb_inner = M_inner*f1
     mb_outer = M_outer*f2
     mb_total = mb_inner + mb_outer

     ! Newtonian acceleration from physical baryonic mass
     gN = mb_total/(r*r)

     y = gN/a0_code

     ! MOND prediction for diagnostic purposes only
     nu = 0.5d0 + sqrt(0.25d0 + 1.d0/y)

     mc = nu*mb_total

     vc      = sqrt(mc/r)
     vc_newt = sqrt(mb_total/r)

     write(99,'(10ES20.10)') &
       r, mb_inner, mb_outer, mb_total, &
       gN, y, nu, mc, vc, vc_newt

  end do

  close(99)

end subroutine print_mond_rotation_table


!========================================================

end module mond_disk_density
