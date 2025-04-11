module declare
  implicit none
  integer, parameter :: j_BounceB = 1, j_elastic = 2, ndir = 18
  integer, parameter :: dirx(0:18) = (/0,+1,-1, 0, 0, 0, 0,+1,-1,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0/)
  integer, parameter :: dirz(0:18) = (/0, 0, 0, 0, 0,+1,-1, 0, 0,+1,-1,+1,-1, 0, 0,-1,+1,-1,+1/)
  integer, parameter :: diry(0:18) = (/0, 0, 0,+1,-1, 0, 0,+1,-1, 0, 0,+1,-1,-1,+1, 0, 0,+1,-1/)
  integer, parameter :: no_slip(0:18) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17/)
  real*8,parameter :: wi(0:18) = (/ 1./3.,1./18,1./18.,1./18.,1./18.,1./18.,1./18., &
                           1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36.,1./36./)
  real*8,parameter :: pi = 4.D0*DATAN(1.D0), rho0 = 1., nu = 1./6., tauinv_Liq = 2./(6.*nu + 1.) 
  real*8,parameter :: omtauinv_Liq = 1. - tauinv_Liq

  real*8,parameter :: Zq_SI = 8.8e6, nu_SI = 1e-6, rho0_SI = 1e3, f0_SI = 5e6, ZL_SI = 716.799 
  real*8,parameter :: delta_f0_SI = sqrt(2 * nu_SI/(2*pi*f0_SI));
  real*8           :: truncation_by_diameter, cell_diam_by_sphere_diam
  real*8           :: f_SI, cell_height_SI, grid_reso_SI, delta_f0,f0
  real*8           :: rad_Sph_SI, rho_Sph_SI, K_Sph_SI, G_Sph_SI, eta_Sph_SI, invtandel_to_Gpp_liq_f0
  real*8           :: x_Sph, y_Sph, z_Sph, rad_Sph, G_Sph, nufac_Sph, K_Sph
  real*8           :: ncycles_mx, delta, omega, omega2by2 
  real*8           :: tauinv_Sph, omtauinv_Sph, tauinv_Liq_m_omega2by2, tauinv_Sph_m_omega2by2
  complex*16 :: omtauinv_Sph_p_Iomega,Iomega
  integer :: i_BounceB_or_Elastic, i_Hard_Soft_Sph, nsteps_mx, ovtorder
  integer :: nx,ny,nz,nLi,nPo,nPoLiq,nPoSph,nPoTop,nPoWalLiq,nPoWalSph,nPoBor,nPoRim,nPo3PL
  
  complex*16, allocatable :: Fx(:), Fy(:), Fz(:)
  complex*16, allocatable :: dr(:),ux(:),uy(:),uz(:),dx(:),dy(:),dz(:),h1(:), h2(:)

  complex*16, allocatable :: uxTop(:), uzTop(:)
  complex*16, allocatable :: Rimux(:),Rimuy(:),Rimuz(:)
  complex*16, allocatable :: FxBor(:,:),FyBor(:,:),FzBor(:,:)
  complex*16, allocatable :: Fx3PLSph(:,:),Fy3PLSph(:,:),Fz3PLSph(:,:)
  complex*16, allocatable :: Fx3PLWal(:,:),Fy3PLWal(:,:),Fz3PLWal(:,:)
  complex*16, allocatable :: FxWalLiq(:,:),FxWalSph(:,:)
  complex*16 :: FxLiqOnWal,FxLiqOnSph,FxSphOnWal,dfc,ZL,Zq,Zliq
  complex*16 dfc_ref,dfc_fin 
    
  integer, allocatable :: Poij(:,:,:),Linj(:,:,:,:),Rimj(:,:,:),Topj(:,:)
  integer, allocatable :: nLiBor(:),nLi3PLSph(:),nLi3PLWal(:)
  integer, allocatable :: PoLiqx(:),PoLiqy(:),PoLiqz(:)
  integer, allocatable :: PoSphx(:),PoSphy(:),PoSphz(:)
  integer, allocatable :: PoTopx(:),PoTopy(:),PoTopz(:)
  integer, allocatable :: PoBorx(:),PoBory(:),PoBorz(:)
  integer, allocatable :: Po3PLx(:),Po3PLy(:),Po3PLz(:)
  integer, allocatable :: PoRimx(:),PoRimy(:),PoRimz(:)
  integer, allocatable :: PoWalLiqx(:),PoWalLiqy(:),PoWalLiqz(:)
  integer, allocatable :: PoWalSphx(:),PoWalSphy(:),PoWalSphz(:)
  integer, allocatable :: permu_Bor(:,:),permu_3PL(:,:)
  integer, parameter :: permu_Top_BBa(1:5 ) = (/        4,      8,        13,12,            18/)
  integer, parameter :: permu_Top_Liq(1:14) = (/0,1,2,3,  5,6,7,  9,10,11,      14,15,16,17   /)
  integer, parameter :: permu_Wal_BBA(1:5 ) = (/      3,      7,       14,      11,      17   /)
  integer, parameter :: permu_Wal_Liq(1:14) = (/0,1,2,  4,5,6,  8,9,10,   12,13,   15,16,   18/)

  integer :: it,nsteps,error
  logical :: B_Save_uxyz
  integer :: beginning,rate,ending
  
  complex*16 :: test,Sauerbrey_Contr
  real*8 :: mass_Sph,target_acc_Hz
end module declare

program FT_LBM
  use declare
  implicit none
  integer :: iovt,x,y,z,iPo,i
  complex * 16 :: dfc_prev 
  logical :: B_accuracy_reached
  
  open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='ux_vs_y.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='uxuy_3D.asc1',     action='write',status='unknown'); close(unit = 20,status='delete');   
  open(unit = 20,file='sim_pars.asc1',    action='write',status='unknown'); close(unit = 20,status='delete');   

  open(unit = 21,file='config2.asc1',      action='read',status='unknown'); 
  read(21,*) i_BounceB_or_Elastic
  print *,'i_BounceB_or_Elastic',   i_BounceB_or_Elastic
  read(21,*) target_acc_Hz
  print *,'target_acc_Hz',   target_acc_Hz
  read(21,*) rad_Sph_SI
  print *,'rad_Sph_SI',   rad_Sph_SI
  read(21,*) cell_height_SI
  print *,'cell_height_SI',   cell_height_SI
  read(21,*) grid_reso_SI
  print *,'grid_reso_SI',   grid_reso_SI
  read(21,*) K_Sph_SI
!  print *,'K_Sph_SI',   K_Sph_SI
  read(21,*) G_Sph_SI 
!  print *,'G_Sph_SI',   G_Sph_SI 
  read(21,*) eta_Sph_SI 
!  print *,'eta_Sph_SI',   eta_Sph_SI 
  read(21,*) rho_Sph_SI
!  print *,'rho_Sph_SI',   rho_Sph_SI
  read(21,*) cell_diam_by_sphere_diam
  print *,'cell_diam_by_sphere_diam',   cell_diam_by_sphere_diam
  read(21,*) truncation_by_diameter
  print *,'truncation_by_diameter',   truncation_by_diameter
  read(21,*) iovt
  print *,'iovt',   iovt
  close(unit = 21);   

  G_Sph_SI = 1e7; eta_Sph_SI = 1e-3; rho_Sph_SI = 1e3
  ncycles_mx = 12; 
  ny = int(cell_height_SI/grid_reso_SI)
  delta_f0 = delta_f0_SI/grid_reso_SI
  f0 = 2.*nu/delta_f0**2 / (2. * pi)
  invtandel_to_Gpp_liq_f0 = G_Sph_SI / (2. * pi * f0_SI * nu_SI * rho0_SI)
  G_Sph = invtandel_to_Gpp_liq_f0 * (2 * pi * f0 * nu * rho0)
  nufac_Sph = (eta_Sph_SI / rho_Sph_SI) / nu_SI
  tauinv_Sph = 2./(6.*nu*nufac_Sph + 1.)
  omtauinv_Sph = 1. - tauinv_sph

  rad_Sph = -1; nx = 1; nz = 1; nPoTop = nx * nz; 

  call Count_and_Allocate
  call make_Lists
  call Set_om_dependent_variables(iovt); 
  nsteps = nsteps_mx
  call init_u0; 
  do it = 1,nsteps_mx
    call compute_VelTop; 
    call Stream_and_Collide; 
    call calc_dfc; 
    if(mod(it,int(nsteps/5)) == 1) then 
      print *,'Jan26,Liq',it,nsteps,real(dfc),aimag(dfc)
    endif  
!    if ( (real(dfc)+aimag(dfc).lt.target_acc_Hz) .and. (errorprev.gt.target_acc_Hz ) ) then 
!      nsteps = it; errorprev = real(dfc)+aimag(dfc)
!    endif  
  enddo
  dfc_ref = dfc;
  B_Save_uxyz = .true. ; call compute_VelTop; call Stream_and_Collide; B_Save_uxyz = .false. ; 
  open(unit = 20,file='ux_vs_y.asc1', action='write',status='unknown', position='append'); 
  do y = 1,ny 
    write (20,*) y*grid_reso_SI*1e9,real(ux(Poij(1,y,1))), aimag(ux(Poij(1,y,1))); 
  enddo  
  close(unit = 20)
  call DeAllocateAll
  rad_Sph = rad_Sph_SI * delta_f0 / delta_f0_SI; 
  nx = int(2. * rad_Sph * cell_diam_by_sphere_diam); 
  if (mod(nx,2).eq.0) nx = nx + 1;
  nz = nx; x_Sph = nx/2.; z_Sph = nz/2.; y_Sph = rad_Sph * (1 - 2. * truncation_by_diameter)
  call system_clock(beginning,rate)
  call Count_and_Allocate
  call make_Lists
  B_Save_uxyz = (i_BounceB_or_Elastic.eq.j_elastic)
  call Set_om_dependent_variables(iovt); 
  call init_u0; 
!  target_acc_Hz = 1e-2
  it = 0  
  B_accuracy_reached = .false.
  do while (it < nsteps_mx .and. .not.B_accuracy_reached )
    it = it + 1  
    call compute_VelTop; 
    if(i_BounceB_or_Elastic.eq.j_BounceB) call sphere_dynamics
    call Stream_and_Collide; 
    call calc_dfc;
    if ( abs(dfc_prev - dfc) < target_acc_Hz ) then 
      print *,'limit reached',it,dfc_prev,dfc
      B_accuracy_reached = .true.
    endif  
    dfc_prev = dfc
    if (nsteps.gt.100) then 
      if(mod(it,int(nsteps/10)) == 1) then 
        print *,it,nsteps,real(dfc-dfc_ref),aimag(dfc-dfc_ref)
      endif 
      if (mod(it,int(nsteps/100)) == 1) then 
        open(unit = 20,file='freq_ring_in.asc1',action='write',status='unknown', position='append'); 
        write (20,*)   ovtorder,it,nsteps,real(dfc-dfc_ref),aimag(dfc-dfc_ref); 
        close(unit = 20);   
      endif  
    endif  
  enddo
  dfc_fin = dfc - dfc_ref
  B_Save_uxyz = .true.; call compute_VelTop; call Stream_and_Collide; B_Save_uxyz = .false.; 
  call system_clock(ending); print *,"Elapsed time is in seconds --->",real(ending - beginning)/real(rate)
  open(unit = 20,file='uxuy_3D.asc1',action='write')
  write(20,"(I5,I5,I5)")nx,ny,nz
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        write(20,"(I5,I5,I5,F12.8,F12.8,F12.8,F12.8,F12.8,F12.8)") x,y,z, &
                                   real(ux(Poij(x,y,z))),aimag(ux(Poij(x,y,z))),&
                                   real(uy(Poij(x,y,z))),aimag(uy(Poij(x,y,z))),&
                                   real(uz(Poij(x,y,z))),aimag(uz(Poij(x,y,z)))
      enddo
    enddo
  enddo
  close(unit = 20)
  print *,dfc_fin/ovtorder
  open(unit = 20,file='sim_pars.asc1',action='write'); 
  write(20,"(I5,I5,I8,I6)") nx,ny,nPo,nsteps
  write(20,*) real(ending - beginning)/real(rate)
  write(20,"(E17.8,E17.8,E17.8)") real((dfc_fin)/ovtorder),aimag((dfc_fin)/ovtorder),grid_reso_SI
  close(unit = 20);   
  call DeAllocateAll
end program FT_LBM

subroutine init_u0
  use declare; implicit none
  integer :: x,y,z,iPo,iLi
  do iPo = 1,nPo
    dr(iPo) = 0
    ux(iPo) = 0
    uy(iPo) = 0
    uz(iPo) = 0
  enddo  
  if (i_BounceB_or_Elastic .eq. j_BounceB) then
    do iPo = 1,nPoSph
      x = PoSphx(iPo)
      y = PoSphy(iPo)
      z = PoSphz(iPo)
      ux(Poij(x,y,z)) = 1
    enddo  
  endif  
  do iLi = 1,nLi
    h1(iLi) = 0
  enddo
end subroutine init_u0

logical function in_Sph(x,y,z) 
  use declare; implicit none
  integer :: x,y,z
  in_Sph = ( (y.ge.1) .and. ((x-x_Sph)**2 + (y-y_Sph)**2 + (z-z_Sph)**2 .lt. rad_Sph**2) )
end function in_Sph

logical function in_Liq(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_mds_not_in_Sph
  B_mds_not_in_Sph = .true.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1)) then 
      B_mds_not_in_Sph = .false.
    endif  
  enddo  
  in_Liq = ((.not. in_Sph(x,y,z)) .and. (y.gt.1) .and. (B_mds_not_in_Sph))
end function in_Liq

logical function on_WalLiq(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_mds_not_in_Sph
  B_mds_not_in_Sph = .true.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1)) then 
      B_mds_not_in_Sph = .false.
    endif  
  enddo  
  on_WalLiq = ((.not. in_Sph(x,y,z)) .and. (y.eq.1) .and. (B_mds_not_in_Sph))
end function on_WalLiq

logical function on_Bor(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_Bor
  B_Bor = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (.not. in_Sph(x,y,z)) .and. in_Sph(xmd,ymd,zmd) .and. (y.gt.1) ) then 
      B_Bor = .true.
    endif  
  enddo  
  on_Bor = B_Bor
end function on_Bor

logical function on_3PL(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_3PL
  B_3PL = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ( (.not. in_Sph(x,y,z)) .and. in_Sph(xmd,ymd,zmd) .and. (ymd.ge.1) .and. (y.eq.1) ) then 
      B_3PL = .true.
    endif  
  enddo  
  on_3PL = B_3PL
end function on_3PL

logical function on_Rim(x,y,z) 
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i
  logical in_Sph,B_Rim
  B_Rim = .false.
  do i = 0,ndir
    xmd = mod(nx+x-dirx(i)-1,nx)+1
    ymd = y-diry(i)
    zmd = mod(nz+z-dirz(i)-1,nz)+1
    if ((.not. in_Sph(xmd,ymd,ymd)) .and. in_Sph(x,y,z)) B_Rim = .true.
  enddo  
  on_Rim = B_Rim
end function on_Rim

subroutine Set_om_dependent_variables(iovt)
  use declare; implicit none
  integer :: iovt
  ovtorder = 2 * iovt - 1; f_SI = ovtorder* f0_SI; 
  delta = delta_f0 / sqrt(float(ovtorder)); 
  omega = 2.*nu/delta**2.; 
  omega2by2 = omega**2/2; Zliq = -conjg(sqrt(rho0*complex(0,1)*omega*rho0*nu))
  tauinv_Liq_m_omega2by2 = tauinv_Liq - omega2by2
  tauinv_Sph_m_omega2by2 = tauinv_Sph - omega2by2
  omtauinv_Sph_p_Iomega = omtauinv_Sph + complex(0.,1.)*omega
  Iomega = complex(0.,1.)*omega
  K_Sph = nu * rho0 * omega * 1e4;
  nsteps_mx = int(ncycles_mx / omega * (ny / delta)**2)
end subroutine Set_om_dependent_variables

subroutine Count_and_Allocate
  use declare; implicit none
  integer :: x,y,z
  logical :: in_Sph,in_Liq,on_WalLiq,on_Bor,on_3PL,on_Rim
  nPo = nx*ny*nz; nLi = nPo * (ndir + 1); 
  nPoLiq = 0; nPoSph = 0; nPoTop = 0; nPoWalLiq = 0; nPoWalSph = 0; nPoBor = 0; nPo3PL = 0; nPoRim = 0;  
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        if ( (in_Liq(x,y,z)) .and. (y.gt.1) .and. (y.lt.ny) )   nPoLiq = nPoLiq + 1        
        if ( on_WalLiq(x,y,z) )                                 nPoWalLiq = nPoWalLiq + 1
        if ( (y.eq.1) .and. in_Sph(x,y,z) )                     nPoWalSph = nPoWalSph + 1
        if (in_Sph(x,y,z))      nPoSph = nPoSph + 1
        if (y.eq.ny)            nPoTop = nPoTop + 1
        if (on_Bor(x,y,z))      nPoBor = nPoBor + 1
        if (on_3PL(x,y,z))      nPo3PL = nPo3PL + 1
        if (on_Rim(x,y,z))      nPoRim = nPoRim + 1
      enddo
    enddo
  enddo
  allocate( h1(nLi),h2(nLi) )
  allocate( dr(nPo),ux(nPo),uy(nPo),uz(nPo) )
  if (i_BounceB_or_Elastic.eq.j_elastic) allocate( dx(nPo),dy(nPo),dz(nPo) )
  allocate( uxTop(nPoTop), uzTop(nPoTop) )
  allocate( Rimux(nPoRim),Rimuy(nPoRim),Rimuz(nPoRim) )
  allocate( FxBor(nPoBor,0:ndir),FyBor(nPoBor,0:ndir),FzBor(nPoBor,0:ndir) )
  allocate( Fx3PLSph(nPo3PL,0:ndir),Fy3PLSph(nPo3PL,0:ndir),Fz3PLSph(nPo3PL,0:ndir) )
  allocate( Fx3PLWal(nPo3PL,0:ndir),Fy3PLWal(nPo3PL,0:ndir),Fz3PLWal(nPo3PL,0:ndir) )
  allocate( FxWalLiq(nPoWalLiq,1:5),FxWalSph(nPoWalSph,1:5) )
  allocate( Poij(nx,ny,nz),Rimj(nx,ny,nz),Topj(nx,nz),Linj(nx,ny,nz,0:ndir) )
  allocate( nLiBor(nPoBor),nLi3PLSph(nPo3PL),nLi3PLWal(nPo3PL) )
  allocate( PoLiqx(nPoLiq),PoLiqy(nPoLiq),PoLiqz(nPoLiq) )
  allocate( PoSphx(nPoSph),PoSphy(nPoSph),PoSphz(nPoSph) )
  allocate( PoTopx(nPoTop),PoTopy(nPoTop),PoTopz(nPoTop) )
  allocate( PoBorx(nPoBor),PoBory(nPoBor),PoBorz(nPoBor) )
  allocate( Po3PLx(nPo3PL),Po3PLy(nPo3PL),Po3PLz(nPo3PL) )
  allocate( PoRimx(nPoRim),PoRimy(nPoRim),PoRimz(nPoRim) )
  allocate( PoWalLiqx(nPoWalLiq),PoWalLiqy(nPoWalLiq),PoWalLiqz(nPoWalLiq) )
  allocate( PoWalSphx(nPoWalSph),PoWalSphy(nPoWalSph),PoWalSphz(nPoWalSph) )
  allocate( permu_Bor(ndir+1,nPoBor),permu_3PL(ndir+1,nPo3PL))
end subroutine Count_and_Allocate

subroutine DeAllocateAll
  use declare; implicit none
  deallocate(h1,h2); 
  deallocate(dr,ux,uy,uz);
  if (i_BounceB_or_Elastic.eq.j_elastic) deallocate(dx,dy,dz)
  deallocate( uxTop, uzTop )
  deallocate( Rimux,Rimuy,Rimuz )
  deallocate( FxBor,FyBor,FzBor)
  deallocate( Fx3PLSph,Fy3PLSph,Fz3PLSph )
  deallocate( Fx3PLWal,Fy3PLWal,Fz3PLWal )
  deallocate( FxWalLiq,FxWalSph )
  deallocate( Poij,Rimj,Linj,Topj )
  deallocate( nLiBor,nLi3PLSph,nLi3PLWal )
  deallocate( PoLiqx,PoLiqy,PoLiqz )
  deallocate( PoSphx,PoSphy,PoSphz )
  deallocate( PoTopx,PoTopy,PoTopz )
  deallocate( PoBorx,PoBory,PoBorz )
  deallocate( Po3PLx,Po3PLy,Po3PLz )
  deallocate( PoRimx,PoRimy,PoRimz )
  deallocate( PoWalLiqx,PoWalLiqy,PoWalLiqz )
  deallocate( PoWalSphx,PoWalSphy,PoWalSphz )
  deallocate( permu_Bor,permu_3PL)
end subroutine DeAllocateAll

subroutine make_Lists
  use declare; implicit none
  integer :: x,y,z,xmd,ymd,zmd,i,i_next
  logical :: in_Sph,in_Liq,on_WalLiq,on_Bor,on_3PL,on_Rim
  nPo = 0; nLi = 0; nPoLiq = 0; nPoSph = 0; nPoTop = 0; nPoWalLiq = 0; nPoWalSph = 0; 
  nPoBor = 0; nPo3PL = 0; nPoRim = 0;  
  do x = 1,nx
    do y = 1,ny
      do z = 1,nz
        nPo = nPo + 1; 
        Poij(x,y,z) = nPo
        do i = 0,ndir
          nLi = nLi + 1 
          Linj(x,y,z,i) = nLi
        enddo

        if ( (in_Liq(x,y,z)) .and. (y.gt.1) .and. (y.lt.ny) ) then
          nPoLiq = nPoLiq + 1        
          PoLiqx(nPoLiq) = x
          PoLiqy(nPoLiq) = y
          PoLiqz(nPoLiq) = z
        endif   

        if ( on_WalLiq(x,y,z) ) then 
          nPoWalLiq = nPoWalLiq + 1
          PoWalLiqx(nPoWalLiq) = x
          PoWalLiqy(nPoWalLiq) = y
          PoWalLiqz(nPoWalLiq) = z
        endif  

        if ( (y.eq.1) .and. in_Sph(x,y,z) ) then 
          nPoWalSph = nPoWalSph + 1
          PoWalSphx(nPoWalSph) = x
          PoWalSphy(nPoWalSph) = y
          PoWalSphz(nPoWalSph) = z
        endif  

        if (in_Sph(x,y,z)) then
          nPoSph = nPoSph + 1
          PoSphx(nPoSph) = x
          PoSphy(nPoSph) = y
          PoSphz(nPoSph) = z
        endif  

        if (y.eq.ny) then 
          nPoTop = nPoTop + 1
          Topj(x,z) = nPoTop
          PoTopx(nPoTop) = x
          PoTopy(nPoTop) = y
          PoTopz(nPoTop) = z
        endif  

        if (on_Bor(x,y,z)) then
          nPoBor = nPoBor + 1
          PoBorx(nPoBor) = x
          PoBory(nPoBor) = y
          PoBorz(nPoBor) = z
          nLiBor(nPoBor) = 0  
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if (in_Sph(xmd,ymd,zmd)) then 
              nLiBor(nPoBor) = nLiBor(nPoBor) + 1
              permu_Bor(nLiBor(nPoBor),nPoBor) = i
            endif
          enddo
          i_next = nLiBor(nPoBor)
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ((ymd.ge.1) .and. (.not.in_Sph(xmd,ymd,zmd)) ) then 
              i_next = i_next + 1
              permu_Bor(i_next,nPoBor) = i
            endif
          enddo
        endif  

        if (on_3PL(x,y,z)) then
          nPo3PL = nPo3PL + 1
          Po3PLx(nPo3PL) = x
          Po3PLy(nPo3PL) = y
          Po3PLz(nPo3PL) = z
          nLi3PLSPh(nPo3PL) = 0  
          nLi3PLWal(nPo3PL) = 0  
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ( (in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1) ) then 
              nLi3PLSPh(nPo3PL) = nLi3PLSPh(nPo3PL) + 1
              permu_3PL(nLi3PLSPh(nPo3PL),nPo3PL) = i
            endif
          enddo
          i_next = nLi3PLSPh(nPo3PL)
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if (ymd.lt.1) then 
              i_next = i_next + 1
              nLi3PLWal(nPo3PL) = nLi3PLWal(nPo3PL) + 1
              permu_3PL(i_next,nPo3PL) = i
            endif
          enddo
          do i = 0,ndir
            xmd = mod(nx+x-dirx(i)-1,nx)+1
            ymd = y-diry(i)
            zmd = mod(nz+z-dirz(i)-1,nz)+1
            if ( ( .not. in_Sph(xmd,ymd,zmd)) .and. (ymd.ge.1) ) then 
              i_next = i_next + 1
              permu_3PL(i_next,nPo3PL) = i
            endif
          enddo
        endif  

        if (on_Rim(x,y,z)) then
          nPoRim = nPoRim + 1
          Rimj(x,y,z) = nPoRim
          PoRimx(nPoRim) = x
          PoRimy(nPoRim) = y
          PoRimz(nPoRim) = z
        endif
      enddo
    enddo
  enddo
end subroutine make_Lists

subroutine compute_VelTop
  use declare; implicit none
  integer :: x,z
  do x = 1,nx
    do z = 1,nz
      uxTop(Topj(x,z)) = -2*( h1(Linj(x,ny,z,7 )) - h1(Linj(x,ny,z,14)) )/(Zliq-1./3.)
      uzTop(Topj(x,z)) = -2*( h1(Linj(x,ny,z,11)) - h1(Linj(x,ny,z,17)) )/(Zliq-1./3.)  ! signs?
    enddo  
  enddo
end subroutine compute_VelTop

subroutine relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: cidotu,heq,drloc,uxloc,uyloc,uzloc,ht(0:ndir)
  cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
  heq = wi(i) * (drloc + 3.*cidotu)
  h2(Linj(x,y,z,i)) = tauinv_Liq_m_omega2by2 * heq + Iomega * ht(i)   !omtauinv = 0
end subroutine relax_Liq

subroutine relax_Sph(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: cidotu,heq,drloc,uxloc,uyloc,uzloc,ht(0:ndir)
  cidotu = dirx(i)*uxloc + diry(i)*uyloc + dirz(i)*uzloc
  heq = wi(i) * (drloc + 3.*cidotu)
  h2(Linj(x,y,z,i)) = omtauinv_Sph_p_Iomega * ht(i) + tauinv_Sph_m_omega2by2 * heq 
end subroutine relax_Sph

subroutine calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
  use declare; implicit none
  integer :: x,y,z,i
  complex*16 :: drloc,ht(0:ndir),uxloc,uyloc,uzloc
  drloc = 0; uxloc = 0; uyloc = 0; uzloc = 0; 
  do i = 0,ndir
    drloc = drloc + ht(i)
    uxloc = uxloc + ht(i) * dirx(i)
    uyloc = uyloc + ht(i) * diry(i)
    uzloc = uzloc + ht(i) * dirz(i)
  enddo  
  if (B_Save_uxyz) then
    dr(Poij(x,y,z)) = drloc
    ux(Poij(x,y,z)) = uxloc
    uy(Poij(x,y,z)) = uyloc
    uz(Poij(x,y,z)) = uzloc
  endif
end subroutine calc_dr_uxyz

subroutine Apply_Elastic_Forces(x,y,z)
  use declare; implicit none
  integer :: x,y,z,xm1,ym1,zm1,xp1,yp1,zp1,i
  complex*16 :: dxWall,Si
  complex*16 :: d_xx_x, d_yy_x, d_zz_x, d_xy_y, d_xz_z
  complex*16 :: d_xx_y, d_yy_y, d_zz_y, d_xy_x, d_yz_z
  complex*16 :: d_xx_z, d_yy_z, d_zz_z, d_xz_x, d_yz_y
  dxWall = 1/(complex(0,1)*omega)
  dx(Poij(x,y,z)) = ux(Poij(x,y,z)) / (complex(0,1)*omega)     
  dy(Poij(x,y,z)) = uy(Poij(x,y,z)) / (complex(0,1)*omega)     
  dz(Poij(x,y,z)) = uz(Poij(x,y,z)) / (complex(0,1)*omega)
  xm1 = mod(nx+x-1-1,nx)+1
  xp1 = mod(nx+x+1-1,nx)+1
  zm1 = mod(nz+z-1-1,nz)+1
  zp1 = mod(nz+z+1-1,nz)+1
  ym1 = y-1
  yp1 = y+1

  d_xx_x = dx(Poij(xp1,y  ,z  )) - 2*dx(Poij(x,y,z)) + dx(Poij(xm1,y  ,z  ))
  d_zz_x = dx(Poij(x  ,y  ,zp1)) - 2*dx(Poij(x,y,z)) + dx(Poij(x  ,y  ,zm1))
  d_xx_y = dy(Poij(xp1,y  ,z  )) - 2*dy(Poij(x,y,z)) + dy(Poij(xm1,y  ,z  ))
  d_zz_y = dy(Poij(x  ,y  ,zp1)) - 2*dy(Poij(x,y,z)) + dy(Poij(x  ,y  ,zm1))
  d_xx_z = dz(Poij(xp1,y  ,z  )) - 2*dz(Poij(x,y,z)) + dz(Poij(xm1,y  ,z  ))
  d_zz_z = dz(Poij(x  ,y  ,zp1)) - 2*dz(Poij(x,y,z)) + dz(Poij(x  ,y  ,zm1))
  d_xz_z = dz(Poij(xp1,y  ,zp1)) - dz(Poij(xp1,y  ,zm1)) - dz(Poij(xm1,y  ,zp1)) + dz(Poij(xm1,y  ,zm1))
  d_xz_x = dx(Poij(xp1,y  ,zp1)) - dx(Poij(xp1,y  ,zm1)) - dx(Poij(xm1,y  ,zp1)) + dx(Poij(xm1,y  ,zm1))

  if (y.gt.1) then 
    d_yy_x = dx(Poij(x  ,yp1,z  )) - 2*dx(Poij(x,y,z)) + dx(Poij(x  ,ym1,z  ))
    d_yy_y = dy(Poij(x  ,yp1,z  )) - 2*dy(Poij(x,y,z)) + dy(Poij(x  ,ym1,z  ))
    d_yy_z = dz(Poij(x  ,yp1,z  )) - 2*dz(Poij(x,y,z)) + dz(Poij(x  ,ym1,z  ))
    d_xy_x = dx(Poij(xp1,yp1,z  )) - dx(Poij(xp1,ym1,z  )) - dx(Poij(xm1,yp1,z  )) + dx(Poij(xm1,ym1,z  ))
    d_xy_y = dy(Poij(xp1,yp1,z  )) - dy(Poij(xp1,ym1,z  )) - dy(Poij(xm1,yp1,z  )) + dy(Poij(xm1,ym1,z  ))
    d_yz_z = dy(Poij(x  ,yp1,zp1)) - dz(Poij(x  ,ym1,zp1)) - dz(Poij(x  ,yp1,zm1)) + dz(Poij(x  ,ym1,zm1))
  endif

  if (y.eq.1) then 
    d_yy_x = dx(Poij(x  ,yp1,z  )) - 2*dx(Poij(x,y,z)) + 0
    d_yy_y = dy(Poij(x  ,yp1,z  )) - 2*dy(Poij(x,y,z)) + 0
    d_yy_z = dz(Poij(x  ,yp1,z  )) - 2*dz(Poij(x,y,z)) + 0
    d_xy_x = dx(Poij(xp1,yp1,z  )) - 0 - dx(Poij(xm1,yp1,z  )) + 0
    d_xy_y = dy(Poij(xp1,yp1,z  )) - 0 - dy(Poij(xm1,yp1,z  )) + 0
    d_yz_z = dy(Poij(x  ,yp1,zp1)) - 0 - dz(Poij(x  ,yp1,zm1)) + 0
  endif  
        
  Fx(Poij(x,y,z)) = (3 * K_Sph - 2 * G_Sph) * (d_xx_x + d_xy_y + d_xz_z) + & 
                    G_Sph * (d_xx_x + d_xy_y + d_xz_z + d_xx_x + d_yy_x + d_zz_x)
  Fy(Poij(x,y,z)) = (3 * K_Sph - 2 * G_Sph) * (d_xy_x + d_yy_y + d_yz_z) + & 
                    G_Sph * (d_xy_x + d_yy_y + d_yz_z + d_yy_y + d_yy_y + d_zz_y)
  Fz(Poij(x,y,z)) = (3 * K_Sph - 2 * G_Sph) * (d_xz_x + d_yz_y + d_zz_z) + &
                    G_Sph * (d_xz_x + d_yz_y + d_zz_z + d_xx_z + d_yy_z + d_zz_z)
  do i = 0,ndir
    Si = (1.-tauinv_Sph/2.)*3.*(dirx(i)*Fx(Poij(x,y,z)) + diry(i)*Fy(Poij(x,y,z)) + dirz(i)*Fz(Poij(x,y,z)))
    h1(Linj(x,y,z,i)) = h1(Linj(x,y,z,i)) + Si
  enddo
end subroutine Apply_Elastic_Forces

subroutine Stream_and_Collide
  use omp_lib
  use declare; implicit none
  integer :: tnr
  integer :: x,y,z,i,iloc,ibar,xmd,ymd,zmd,iPo,iLi,iRim
  complex*16 :: ht(0:ndir),drloc,uxloc,uyloc,uzloc,cidotu
  do iPo = 1,nPoBor
    do iloc = 0,ndir
      FxBor(iPo,iloc) = 0
      FyBor(iPo,iloc) = 0
      FzBor(iPo,iloc) = 0
    enddo
  enddo  

  do iPo = 1,nPoWalLiq
    x = PoWalLiqx(iPo)
    y = PoWalLiqy(iPo)
    z = PoWalLiqz(iPo)
    do iloc = 1,5
      i = permu_Wal_BBa(iloc)
      ibar = no_slip(i)
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * dirx(ibar)
      FxWalLiq(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
    enddo  
    do iloc = 1,14
      i = permu_Wal_Liq(iloc)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  
!  !$omp parallel private(iPo,i,x,y,z,drloc,uxloc,uyloc,xmd,ymd,zmd,uzloc,ht)
!    !$omp do
      do iPo = 1,nPoLiq
        x = PoLiqx(iPo)
        y = PoLiqy(iPo)
        z = PoLiqz(iPo)
        do i = 0,ndir
          xmd = mod(nx+x-dirx(i)-1,nx)+1
          ymd = y-diry(i)
          zmd = mod(nz+z-dirz(i)-1,nz)+1
          ht(i) = h1(Linj(xmd,ymd,zmd,i))
        enddo
        call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
        do i = 0,ndir
          call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
        enddo  
      enddo  
!    !$omp end do
!  !$omp end parallel
  if (i_BounceB_or_Elastic.eq.j_elastic) then
    do iPo = 1,nPoWalSph
      x = PoWalSphx(iPo)
      y = PoWalSphy(iPo)
      z = PoWalSphz(iPo)
      do iloc = 1,5
        i = permu_Wal_BBa(iloc)
        ibar = no_slip(i)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * dirx(ibar)
        FxWalSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
      enddo  
      do iloc = 1,14
        i = permu_Wal_Liq(iloc)
        xmd = mod(nx+x-dirx(i)-1,nx)+1
        ymd = y-diry(i)
        zmd = mod(nz+z-dirz(i)-1,nz)+1
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      enddo  
      call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
      do i = 0,ndir
        call relax_Sph(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
      enddo  
      call Apply_Elastic_Forces(x,y,z)
    enddo
  endif  

  if (i_BounceB_or_Elastic.eq.j_elastic) then
    do iPo = 1,nPoSph
      x = PoSphx(iPo)
      y = PoSphy(iPo)
      z = PoSphz(iPo)
      do i = 0,ndir
        xmd = mod(nx+x-dirx(i)-1,nx)+1
        ymd = y-diry(i)
        zmd = mod(nz+z-dirz(i)-1,nz)+1
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      enddo  
      call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
      do i = 0,ndir
        call relax_Sph(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
      enddo  
      call Apply_Elastic_Forces(x,y,z)
    enddo  
  endif  

  do iPo = 1,nPoBor
    x = PoBorx(iPo)
    y = PoBory(iPo)
    z = PoBorz(iPo)
    do iloc = 1,nLiBor(iPo)
      i = permu_Bor(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      if (i_BounceB_or_Elastic.eq.j_BounceB) then
        ibar = no_slip(i)
        iRim = Rimj(xmd,ymd,zmd)
        cidotu = dirx(ibar)*Rimux(iRim) + diry(ibar)*Rimuy(iRim) + diry(ibar)*Rimuz(iRim)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
        FxBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
        FyBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * diry(i)
        FzBor(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirz(i)
      else   
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      endif
    enddo  
    do iloc = nLiBor(iPo)+1,ndir+1
      i = permu_Bor(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  

  do iPo = 1,nPo3PL
    x = Po3PLx(iPo)
    y = Po3PLy(iPo)
    z = Po3PLz(iPo)
    do iloc = 1,nLi3PLSPh(iPo)
      i = permu_3PL(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      if (i_BounceB_or_Elastic.eq.j_BounceB) then
        ibar = no_slip(i)
        iRim = Rimj(xmd,ymd,zmd)
        cidotu = dirx(ibar)*Rimux(iRim) + diry(ibar)*Rimuy(iRim) + diry(ibar)*Rimuz(iRim)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
        Fx3PLSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
        Fy3PLSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * diry(i)
        Fz3PLSph(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirz(i)
      else   
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      endif
    enddo  
    do iloc = nLi3PLSPh(iPo) + 1, nLi3PLSPh(iPo) + nLi3PLWal(iPo)
      i = permu_3PL(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      if (i_BounceB_or_Elastic.eq.j_BounceB) then
        ibar = no_slip(i)
        ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * dirx(ibar)
        Fx3PLWal(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirx(i)
        Fy3PLWal(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * diry(i)
        Fz3PLWal(iPo,iloc) = (ht(i) + h1(Linj(x,y,z,ibar))) * dirz(i)
      else   
        ht(i) = h1(Linj(xmd,ymd,zmd,i))
      endif
    enddo  
    do iloc = nLi3PLSPh(iPo) + nLi3PLWal(iPo) + 1,ndir+1
      i = permu_3PL(iloc,iPo)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo  

  do iPo = 1,nPoTop
    x = PoTopx(iPo)
    y = PoTopy(iPo)
    z = PoLiqz(iPo)
    do iloc = 1,5
      i = permu_Top_BBa(iloc)
      ibar = no_slip(i)
      cidotu = dirx(ibar)*uxTop(Topj(x,z)) + dirz(ibar)*uzTop(Topj(x,z))
      ht(i) = h1(Linj(x,y,z,ibar)) - 6 * wi(ibar) * cidotu
    enddo  
    do iloc = 1,14
      i = permu_Top_Liq(iloc)
      xmd = mod(nx+x-dirx(i)-1,nx)+1
      ymd = y-diry(i)
      zmd = mod(nz+z-dirz(i)-1,nz)+1
      ht(i) = h1(Linj(xmd,ymd,zmd,i))
    enddo  
    call calc_dr_uxyz(x,y,z,ht,drloc,uxloc,uyloc,uzloc)
    do i = 0,ndir
      call relax_Liq(x,y,z,i,drloc,uxloc,uyloc,uzloc,ht)
    enddo  
  enddo

  do iLi = 1,nLi
    h1(iLi) = h2(iLi)
  enddo
end subroutine Stream_and_Collide

subroutine sphere_dynamics
  use declare; implicit none
  integer :: iPo,iloc
  real :: al 
  al = truncation_by_diameter  ! ???
  mass_Sph = (rho_Sph_SI/rho0_SI) * rho0 * rad_Sph**3. * 4.*pi/3 *(al-1.)**2 * (1 + 2*al)
  !4/3 (-1 + al)^2 (1 + 2 al) \[Pi]
  
  FxLiqOnSph = 0.
  do iPo = 1,nPoBor
    do iloc = 1,nLiBor(iPo)
      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
    enddo
  enddo
  do iPo = 1,nPo3PL
    do iloc = 1,nLi3PLSph(iPo)
      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
    enddo
  enddo
  FxSphOnWal = FxLiqOnSph - complex(0.,1.) * omega * mass_Sph 
  do iPo = 1,nPoRim
    Rimux(iPo) = 1;
    Rimuy(iPo) = 0;
    Rimuz(iPo) = 0;
  enddo  
end subroutine sphere_dynamics

subroutine calc_dfc
  use declare; implicit none
  integer :: iPo,iloc
  FxLiqOnWal = 0.
  do iPo = 1,nPoWalLiq
    do iloc = 1,5
      FxLiqOnWal = FxLiqOnWal + FxWalLiq(iPo,iloc)     
    enddo
  enddo
    do iPo = 1,nPo3PL
    do iloc = nLi3PLSph(iPo) + 1, nLi3PLSph(iPo) + nLi3PLWal(iPo)
      FxLiqOnWal = FxLiqOnWal + Fx3PLWal(iPo,iloc)
    enddo
  enddo

  if (i_BounceB_or_Elastic.eq.j_elastic) then
    do iPo = 1,nPoWalSph
      do iloc = 1,5
        FxLiqOnWal = FxLiqOnWal + FxWalSph(iPo,iloc)     
      enddo
    enddo  
  else
    FxLiqOnWal = FxLiqOnWal + FxSphOnWal
  endif  
  ZL = (FxLiqOnWal)/(nx*nz)
  Zq = (rho0/rho0_SI)*(nu/nu_SI)*(delta_f0_SI/delta_f0)*Zq_SI
  Sauerbrey_Contr = complex(0,1) * FxSphOnWal/(nx*nz) * f0_SI /(pi*Zq)
  Sauerbrey_Contr = -conjg(Sauerbrey_Contr)
  dfc = f0_SI * complex(0,1)/(pi*Zq) * ZL
  dfc = -conjg(dfc)
end subroutine calc_dfc

!subroutine Update_Top_Adv
!  use declare; implicit none
!  integer :: x,z,iqx,iqz,iqxmax,iqzmax
!  real*8 :: qx,qz,q,phi
!!  complex*16 :: uxrot,uxqzbar,uzqxbar,uzqzbar
!!  complex*16 :: sigxyqxbar,sigxyqzbar,sigzyqxbar,sigzyqzbar,sigyyqxbar,sigyyqzbar  
!  complex*16 :: uxq,uyq,sigxyq,sigyyq,kSlow,kFast,sigxy,sigzy,sigyy
!  complex*16, allocatable :: uxSlow(:,:),uxFast(:,:),uzSlow(:,:),uzFast(:,:) 
!  iqxmax = 1; iqxmax = 1; 
!  allocate( uxSlow(-iqxmax:iqxmax,-iqzmax:iqzmax) )
!  allocate( uzSlow(-iqxmax:iqxmax,-iqzmax:iqzmax) )
!  allocate( uxFast(-iqxmax:iqxmax,-iqzmax:iqzmax) )
!  allocate( uzFast(-iqxmax:iqxmax,-iqzmax:iqzmax) )
!  do iqx = -iqxmax,iqxmax 
!    qx = iqx * pi / nx
!    do iqz = -iqzmax,iqzmax 
!      qz = iqz * pi / nx
!      phi = atan(qz,qz) 
!      uxq = 0; uyq = 0; sigxyq = 0; sigyyq = 0; 
!      q = (qx**2 + qz**2)**(1/2)
!      do x = 1,nx
!        do z = 1,nz
!          uxq = uxq + ( ux(Poij(x,ny,z)) * cos(phi) + uz(Poij(x,ny,z)) * sin(phi) )  * &
!            exp( complex(0,1) * ( x*cos(phi) + z*sin(phi) ) * q) ; 
!          uyq = uyq + uy(Poij(x,ny,z)) * &
!            exp( complex(0,1) * ( x*cos(phi) + z*sin(phi) ) * q) ; 
!          sigxy = 1;
!          sigzy = 1; 
!          sigyy = 1;
!        enddo
!      enddo  
!    enddo
!  enddo
!end subroutine Update_Top_Adv

!subroutine Update_Sphere_Motion
!  use declare; implicit none
!  complex * 16 :: uSphx,uSphy,uSphz,OmSphx,OmSphy,OmSphz,&
!               FyLiqOnSph,FzLiqOnSph, &  ! FxLiqOnSph is global
!    TxLiqOnSph,TyLiqOnSph,TzLiqOnSph, &
!    Expan,Stres_Expan, &
!    Shear_xy,Shear_xz,Shear_yz,Shear_zmx,Shear_ymxz, &
!    Stres_xy,Stres_xz,Stres_yz,Stres_zmx,Stres_ymxz, &
!    z_CM,a_CR,SpringConst,BendStiffn, &
!    External_Stres(3,3)
!  integer :: i1,i2    
!    
!  FxLiqOnSph = 0.; FyLiqOnSph = 0.; FzLiqOnSph = 0.
!  do i1 = 1,3
!    do i2 = 1,3
!      External_Stres(i1,i2) = 0.
!    enddo
!  enddo   ! pxr ?
!  do iPo = 1,nPoBor
!    do iloc = 1,nLiBor(iPo)
!      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
!    enddo
!  enddo
!  do iPo = 1,nPo3PL
!    do iloc = 1,nLi3PLSph(iPo)
!      FxLiqOnSph = FxLiqOnSph + FxBor(iPo,iloc)
!    enddo
!  enddo
!end subroutine Update_Sphere_Motion
!!var idim : Byte;
!    uBC,
!    ExRC,SxyC,SxzC,SyzC,SxmyC,SzmxyC,
!    FexRC,FSxyC,FSxzC,FSyzC,FSxmyC,FSzmxyC,
!    FLinkexR,FLinkSxy,FLinkSxz,FLinkSyz,FLinkSxmy,FLinkSzmxy,
!    FelasexR,FelasSxy,FelasSxz,FelasSyz,FelasSxmy,FelasSzmxy,
!    FtotexR,FtotSxy,FtotSxz,FtotSyz,FtotSxmy,FtotSzmxy,
!    DelExRC,DelSxyC,DelSxzC,DelSyzC,DelSxmyC,DelSzmxyC,
!    ExRCnew,SxyCnew,SxzCnew,SyzCnew,SxmyCnew,SzmxyCnew,
!    adjfacuC,adjfacRR,adjfacEx : TComplex;
!    adjdeluC,adjmaguC,adjdelRR,adjmagRR,adjdelEx,adjmagEx : Double;
!begin
!  adjdeluC := Pi/4;
!  adjdelRR := Pi/4;
!  adjdelEx := Pi/4;
!  adjmaguC := Main.adj_fac_gl * Main.i_updr;
!  adjmagRR := 0.5 * Main.adj_fac_gl * Main.i_updr;
!  adjmagEx := 0.5 * Main.adj_fac_gl * Main.i_updr;
!  adjfacuC.real := adjmaguC * Cos(adjdeluC); adjfacuC.imag := -adjmaguC * Sin(adjdeluC);
!  adjfacRR.real := adjmagRR * Cos(adjdelRR); adjfacRR.imag := -adjmagRR * Sin(adjdelRR);
!  adjfacEx.real := adjmagEx * Cos(adjdelEx); adjfacEx.imag := -adjmagEx * Sin(adjdelEx);
!  with Compl do begin
!    uBC := CPlx(uB[1],uB[2]);
!    velB[1] := CPlx(uB[1],uB[2]);               velB[2] := C0;                              velB[3] := C0;
!    velC[1] := CPlx(uC[1,i_sph],uC[2,i_sph]);   velC[2] := CPlx(vC[ 1,i_sph],vC[ 2,i_sph]); velC[3] := CPlx( wC[1,i_sph], wC[2,i_sph]);
!    RRC[ 1] := CPlx(RRx[1,i_sph],RRx[2,i_sph]); RRC[2]  := CPlx(RRy[1,i_sph],RRy[2,i_sph]); RRC[3]  := CPlx(RRz[1,i_sph],RRz[2,i_sph]);
!    FC[1]   := CPlx(Fx[1,i_sph],Fx[2,i_sph]);   FC[2]   := CPlx(Fy[ 1,i_sph],Fy[ 2,i_sph]); FC[3]   := CPlx( Fz[1,i_sph], Fz[2,i_sph]);
!    TC[1]   := CPlx(Tx[1,i_sph],Tx[2,i_sph]);   TC[2]   := CPlx(Ty[ 1,i_sph],Ty[ 2,i_sph]); TC[3]   := CPlx( Tz[1,i_sph], Tz[2,i_sph]);
!
!
!    ExRC   := CPlx(ExR[1,i_sph],ExR[2,i_sph]);
!    ExRC   := C0;
!    SxyC   := CPlx(Sxy[1,i_sph],Sxy[2,i_sph]);
!    SxzC   := CPlx(Sxz[1,i_sph],Sxz[2,i_sph]);
!    SyzC   := CPlx(Syz[1,i_sph],Syz[2,i_sph]);
!    SxmyC  := CPlx(Sxmy[1,i_sph],Sxmy[2,i_sph]);
!    SzmxyC := CPlx(Szmxy[1,i_sph],Szmxy[2,i_sph]);
!
!    FexRC   := CPlx(FexR[1,i_sph],FexR[2,i_sph]);
!    FSxyC   := CPlx(FSxy[1,i_sph],FSxy[2,i_sph]);
!    FSxzC   := CPlx(FSxz[1,i_sph],FSxz[2,i_sph]);
!    FSyzC   := CPlx(FSyz[1,i_sph],FSyz[2,i_sph]);
!    FSxmyC  := CPlx(FSxmy[1,i_sph],FSxmy[2,i_sph]);
!    FSzmxyC := CPlx(FSzmxy[1,i_sph],FSzmxy[2,i_sph]);
!
!    if (ac > 0) and (Sph_zG < Sph_R) then begin
!      velSB[1].real := velC[1].real - Sph_zCM*(RRC[2].real+SxzC.real);
!      velSB[1].imag := velC[1].imag - Sph_zCM*(RRC[2].imag+SxzC.imag);
!      velSB[2].real := velC[2].real - Sph_zCM*(RRC[1].real+SyzC.real);
!      velSB[2].imag := velC[2].imag - Sph_zCM*(RRC[1].imag+SyzC.imag);
!      velSB[3].real := velC[3].real - Sph_zCM*(exRC.real+SzmxyC.real);;
!      velSB[3].imag := velC[3].imag - Sph_zCM*(exRC.imag+SzmxyC.imag);;
!      for idim  := 1 to 3 do begin
!        velBmS[idim] := CMinus(velB[idim],velSB[idim]);
!        case idim of
!          1,2 : begin
!            FLinkLoc[idim] := CTimes(kapSbyiom,velBmS[idim]);
!            TLinkLoc[idim] := CTimes(kapBbyiom,RRC[idim]);
!          end;
!          3   : begin
!            FLinkLoc[idim] := CTimes(kapNbyiom,velBmS[idim]);
!            TLinkLoc[idim] := CTimes(kapTbyiom,RRC[idim]);
!          end;
!        end;
!      end;
!      FLinkexR := CTimes(kapNbyiom,CDTimes(-Sph_zCM,velBmS[3]));
!      FLinkexR := C0;
!      FLinkSxy := C0;
!      FLinkSxz := CTimes(kapSbyiom,CDTimes(-Sph_zCM,velBmS[1]));
!      FLinkSyz := CTimes(kapSbyiom,CDTimes(-Sph_zCM,velBmS[2]));
!      FLinkSxmy := C0;
!      FLinkSzmxy := CTimes(kapNbyiom,CDTimes(-Sph_zCM,velBmS[3]));
!
!      FLink[1] := CPlus(FLinkLoc[1],CDTimes(1/Sph_zCM,TLinkLoc[2]));
!      FLink[2] := CPlus(FLinkLoc[2],CDTimes(1/Sph_zCM,TLinkLoc[1]));
!      FLink[3] :=       FLinkLoc[3];
!      TLink[1] := CPlus(TLinkLoc[1],CDTimes(Sph_zCM,FLinkLoc[2]));
!      TLink[2] := CPlus(TLinkLoc[2],CDTimes(Sph_zCM,FLinkLoc[1]));
!      TLink[3] :=       TLinkLoc[3];
!    end else begin
!      for idim  := 1 to 3 do begin FLink[idim] := C0; TLink[idim] := C0; end;
!      FLinkexR := C0; FLinkSxy := C0; FLinkSxz := C0; FLinkSyz := C0; FLinkSxmy := C0; FLinkSzmxy := C0;
!    end;
!    FelasexR   := CTimes(CDTimes(1,KtVbyiom),exRC);
!    FelasSxy   := CTimes(CDTimes(1,GtVbyiom),SxyC);
!    FelasSxz   := CTimes(CDTimes(1,GtVbyiom),SxzC);
!    FelasSyz   := CTimes(CDTimes(1,GtVbyiom),SyzC);
!    FelasSxmy  := CTimes(CDTimes(1,GtVbyiom),SxmyC);
!    FelasSzmxy := CTimes(CDTimes(1,GtVbyiom),SzmxyC);
!
!    for idim  := 1 to 3 do begin
!      Ftot[idim] := CPlus(FC[idim],FLink[idim]);
!      Ttot[idim] := CPlus(TC[idim],TLink[idim]);
!    end;
!    FtotexR := CPlus(CPlus(FexRC,FLinkexR),FelasexR);
!    FtotSxy := CPlus(CPlus(FSxyC,FLinkSxy),FelasSxy);
!    FtotSxz := CPlus(CPlus(FSxzC,FLinkSxz),FelasSxz);
!    FtotSyz := CPlus(CPlus(FSyzC,FLinkSyz),FelasSyz);
!    FtotSxmy := CPlus(CPlus(FSxmyC,FLinkSxmy),FelasSxmy);
!    FtotSzmxy := CPlus(CPlus(FSzmxyC,FLinkSzmxy),FelasSzmxy);
!
!    for idim  := 1 to 3 do DelvelC[idim] := CMinus(CDiv(Ftot[idim],iomMass),velC[idim]);
!    DelRRC[1] := CMinus(CDiv(Ttot[1],iomIxx),RRC[1]);
!    DelRRC[2] := CMinus(CDiv(Ttot[2],iomIyy),RRC[2]);
!    DelRRC[3] := CMinus(CDiv(Ttot[3],iomIzz),RRC[3]);
!
!    DelExRC   := CMinus(CDiv(FtotexR,iomIall),ExRC);
!    DelSxyC   := CMinus(CDiv(FtotSxy,iomIzz ),SxyC);
!    DelSxzC   := CMinus(CDiv(FtotSxz,iomIyy ),SxzC);
!    DelSyzC   := CMinus(CDiv(FtotSyz,iomIzz ),SyzC);
!    DelSxmyC  := CMinus(CDiv(FtotSxmy,iomIzz),SxmyC);
!    DelSzmxyC := CMinus(CDiv(FtotSzmxy,iomIzmxy),SzmxyC);
!
!    for idim  := 1 to 3 do begin
!      velCnew[idim] := CPlus(velC[idim],CTimes(adjfacuC,DelvelC[idim]));
!      RRCnew[idim]  := CPlus(RRC[idim],CTimes(adjfacRR,DelRRC[idim]));
!    end;
!
!    ExRCnew   := CPlus(ExRC,CTimes(adjfacEx,DelExRC));
!    ExRCnew   := C0;
!    SxyCnew   := CPlus(SxyC,CTimes(adjfacRR,DelSxyC));
!    SxzCnew   := CPlus(SxzC,CTimes(adjfacRR,DelSxzC));
!    SyzCnew   := CPlus(SyzC,CTimes(adjfacRR,DelSyzC));
!    SxmyCnew  := CPlus(SxmyC,CTimes(adjfacRR,DelSxmyC));
!    SzmxyCnew := CPlus(SzmxyC,CTimes(adjfacRR,DelSzmxyC));
!
!    uC[1,i_sph] := velCnew[1].real; uC[2,i_sph] := velCnew[1].imag;
!    vC[1,i_sph] := velCnew[2].real; vC[2,i_sph] := velCnew[2].imag;
!    wC[1,i_sph] := velCnew[3].real; wC[2,i_sph] := velCnew[3].imag;
!    RRx[1,i_sph] := RRCnew[1].real; RRx[2,i_sph] := RRCnew[1].imag;
!    RRy[1,i_sph] := RRCnew[2].real; RRy[2,i_sph] := RRCnew[2].imag;
!    RRz[1,i_sph] := RRCnew[3].real; RRz[2,i_sph] := RRCnew[3].imag;
!
!    ExR[1,i_sph] := ExRCnew.real; ExR[2,i_sph] := ExRCnew.imag;
!    Sxy[1,i_sph] := SxyCnew.real; Sxy[2,i_sph] := SxyCnew.imag;
!    Sxz[1,i_sph] := SxzCnew.real; Sxz[2,i_sph] := SxzCnew.imag;
!    Syz[1,i_sph] := SyzCnew.real; Syz[2,i_sph] := SyzCnew.imag;
!    Sxmy[1,i_sph] := SxmyCnew.real; Sxmy[2,i_sph] := SxmyCnew.imag;
!    Szmxy[1,i_sph] := SzmxyCnew.real; Szmxy[2,i_sph] := SzmxyCnew.imag;
!    Lists.Update_vel_Spheres;
!
