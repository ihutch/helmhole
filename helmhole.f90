! Construct electron hole equilibria on an rz mesh using a
! Modified Helmholtz (flat-trap shielding) term plus 
! a trapped particle deficit represented by the offset power model. 
! Uses the fishpack library subroutine sepeli to solve the Helmholtz equation. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module helmhole
  integer, parameter :: nz=100,nr=100,nclv=24,nprof=10
  integer, parameter :: nwork=nz*nr
  real, dimension(nwork) :: WORK,pwork
  real, dimension(0:nr) :: r,BDA,BDB
  real, dimension(0:nz) :: z,BDC,BDD
  real, dimension(0:nr,0:nz) :: phi,rhominus,USOL,rhominusnew,Ez,Er
  real, dimension(nclv) :: zclv
  real, dimension(0:nz) :: pfit
  real, dimension(0:nz,nprof) :: phinprof,rhonprof
  real, dimension(0:nr) :: Erbypsi,Ezp,Erp,Erzr
  real phiguarded(-1:nz+1,-1:nr+1)
! Scalar parameter values
  real :: zmin=0.,zmax=10.,rmin=0.,rmax=10.,dz,dr
  real :: dlength=1., um=0.0, Omega=1.
  real :: psimax=.1, wj=-.02, alpha=.5, wjmod, psimaxd
  real :: rt=5. , sf=1., sh=0., vperpfac=1., psiupfac=0.0, cgfac=1.
  integer :: niter=40, ips=0
  logical :: lradial=.true.,lreadf=.false.,lrprof=.true.,lzlog=.false.
  logical :: lsparselabels=.true.,lzprof=.false.
  character*30 arg
  character cyl
  real :: delta
! Transverse Profiles
  real, dimension(0:nr) :: psi0r, psidr, wjr, wjmodr, fpsir, del2perp, Fperp
  real, dimension(0:nr) :: Erbyphi
! Boundary condition specifiers. 3 both mixed| 4 slope/fixed. 
  integer :: IORDER=2,MBDCND=3,NBDCND=4,NBDeli=3,INTL=0
  real :: ASEPX,BSEPX,GSEPY,XNU
contains
  ! Mesh and Boundary condition setting.
  subroutine initrzbd
    do i=0,nr
       r(i)=rmin+(rmax-rmin)*(i)/(nr)
       USOL(i,nz)=0. ! zero at zmax
       USOL(i,0)=0.
       BDC(i)=0.     ! zero slope at zmin
       GSEPY=0.
       BDD(i)=0.     ! mixed sum value at zmax
       XNU=dlength   ! logarithmic derivative.
    enddo
    do j=0,nz
       z(j)=zmin+(zmax-zmin)*(j)/(nz)
       BDA(j)=0. ! Value of mixed bc sum
       ASEPX=0.  ! Coefficient of phi (phi' coef is 1, so phi'=0.)
       BDB(j)=0.
       BSEPX=dlength  ! at B dln(p)/dr=-1
    enddo
    WORK(1)=nwork
  end subroutine initrzbd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine setdlength 
! Take the screening length to be ~ sqrt of slope of refdenapprox
    dp=1.e-3*psimax
    dn=refdenapprox(dp,um)-1.
    dlength=1./sqrt(dn/dp)
    write(*,*)'Um=',um,' dn=',dn,' Dlength=',dlength
  end subroutine setdlength
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine parsearguments
    do i=1,iargc()   ! Parse cmdline arguments.
       call getarg(i,arg)
       if(arg(1:2).eq.'-w')read(arg(3:),*)wj
       if(arg(1:2).eq.'-a')read(arg(3:),*)alpha
       if(arg(1:2).eq.'-u')read(arg(3:),*)um
       if(arg(1:2).eq.'-O')read(arg(3:),*)Omega
       if(arg(1:2).eq.'-f')lreadf=.true.
       if(arg(1:3).eq.'-rp')lrprof=.not.lrprof
       if(arg(1:3).eq.'-ia')read(arg(4:),*)psiupfac
       if(arg(1:3).eq.'-zm')read(arg(4:),*)zmax
       if(arg(1:3).eq.'-zp')lzprof=.not.lzprof
       if(arg(1:3).eq.'-rm')read(arg(4:),*)rmax
       if(arg(1:3).eq.'-rt')read(arg(4:),*)rt
       if(arg(1:3).eq.'-sf')read(arg(4:),*)sf
       if(arg(1:3).eq.'-sh')read(arg(4:),*)sh
       if(arg(1:3).eq.'-ps')then
          read(arg(4:),*)ips
       elseif(arg(1:2).eq.'-p')then
          read(arg(3:),*)psimax
       endif
       if(arg(1:2).eq.'-h')goto 1
       if(arg(1:3).eq.'-vp')read(arg(4:),*)vperpfac
    enddo
    if(lradial)then; cyl='r'; else; cyl='y'; endif
    call setdlength   
    call pfset(ips)
    return
1   write(*,*)'Usage: helmhole [-p<value> -w -a -u -zm -rm -ps -h ...]'
    write(*,*)'Switch     Parameter          [Current]'
    write(*,101)'-p   set peak potential desired [',psimax,']'
    write(*,101)'-w   set negative energy offset [',wj,']'
    write(*,101)'-a   set power alpha            [',alpha,']'
    write(*,101)'-u   set hole speed             [',um,']'
    write(*,101)'-zm  set domain z maximum       [',zmax,']'
    write(*,101)'-rm  set domain r maximum       [',rmax,']'
    write(*,101)'-rt  set radial flat-top length [',rt,']'
    write(*,101)'-sh  set extra r-curvature      [',sh,']'
    write(*,101)'-sf  override shielding length  [',sf,']'
    write(*,103)'-ps  set plot-files type e.g. 3 [',ips,' ]'
    write(*,102)'-rp  toggle parallel functional [',lrprof,']'
    write(*,*)'r-shape =(1.+sg)/(cosh(sf*y(j)/dl)+sg*(1+sh*y(j)**2); sg=exp(rt)'
    write(*,*)'Parameter values are trailing part of switch without whitespace'
    write(*,'(a)')'-h -?  Print this help.'
    call exit
101 format(a,f5.2,a,f5.2)
102 format(a,l5,a,l5)
103 format(a,i4,a,i4)
  end subroutine parsearguments
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine initwbpsi
! Initialize the profiles of wb and psi. r and z already initialized.
    real :: vflatapprox
    integer :: i,j
    external vflatapprox
    dz=(zmax-zmin)/nz                    ! z0 is z=0
    dz2inv=1./dz**2
    dr=(rmax-rmin)/nr
    dr2inv=1./dr**2
! radial-shape parameters of psi    
    sg=exp(rt) ! Argument is approximate transverse radius.
! sg controls flattop, sh is flattop minus on-axis curvature: -d2p/dr2/p
! sf is 1/dlength (by default).
    do j=0,nr
       sy=(1.+sg)/(cosh(sf*r(j)/dlength)+sg*(1.+sh*r(j)**2))  !r-shape of psi
       psi0r(j)=psimax*sy  
       wjr(j)=wj*sy         
    enddo
    do j=1,nr-1  ! delperp^2 cylindrical normalized phi.
       del2perp(j)=((r(j)+r(j+1))*(psi0r(j+1)-psi0r(j)) -  &
            (r(j-1)+r(j))*(psi0r(j)-psi0r(j-1)))/(2.*r(j))
    enddo
    del2perp(0)=4.*(psi0r(1)-psi0r(0))   ! assumes r(0)=0.
    del2perp(nr)=del2perp(nr-1)*psi0r(nr)/psi0r(nr-1)  
    del2perp=del2perp*dr2inv/psi0r  ! (1/psi(y))(d/dy)^2(psi)
! Now fpsir and phi initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Originally my hand notes used units of velocity u=v/sqrt(2T/m) and f(u).
! In the ftfpowermodel write-up however, it is a distraction to introduce
! that differently scaled velocity. So only v/sqrt(T/m) units are used.
! f(v)dv=f(u)du, so f(v)=f(u)(du/dv)=f(u)/sqrt(2). All that is needed
! to convert the whole code to v-units is to multiply the definition of
! G by sqrt(2). That has been done so this code is now v-units not u-units.
! Everything else is a function of phi, so there is no problem with that.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    G=sqrt(3.1415926)*gamma(alpha+1.)/(2.*gamma(alpha+1.5))  !u-units
    G=sqrt(3.1415926/2.)*gamma(alpha+1.)/(gamma(alpha+1.5))  !v-units
    do j=0,nr
       Fperp(j)=-(2*Vflatapprox(psi0r(j),um)/psi0r(j)**2+vperpfac*del2perp(j))
       fpsir(j)=-(alpha+3/2.)*psi0r(j)**2 &
            /(4.*G*(psi0r(j)+wjr(j))**1.5)*Fperp(j)
!       write(*,*)j,Fperp(j),'fpsir=',fpsir(j)
       fpsir(j)=min(fpsir(j),0.)
       do i=0,nz
          USOL(j,i)=psi0r(j)/cosh(z(i)/3.)**2    ! Initial z shape arbitrary.
       enddo
    enddo
    psidr=psi0r ! psidr is psi in rho to get psi0dr.
    psimaxd=psimax
!    write(*,*)sf,sh,rt,sg
    if(.not.lzprof)then
    call charsize(0.018,0.018)
    call autoplot(r,psi0r,nr)
    call axlabels('r','')
    call boxtitle('Prescribed Radial Profiles')
    call jdrwstr(wx2nx(1.),wy2ny(-wjr(0))+.03,'-W!d!Bjd!@!d(r)',1.)
    ipr=int(nr/3.)
    call jdrwstr(wx2nx(r(ipr)),wy2ny(psi0r(ipr)),' !Ay!B!dd!d!@(r)',1.)
    call dashset(2)
    call polyline(r,-wjr,nr)
    call dashset(0)
    call pltend
    endif
! Very rough psi initialization because it only affects delta.
    do j=0,nr
       phi(j,:)=psi0r(j)/(cosh(z/2.)**2)
    enddo
  end subroutine initwbpsi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initntilde
! Initialize the potential as the solution of the interior Helmholtz
! equation expressed as a sum of harmonic solutions, out to the value
! -Wj which defines the join contour. ntilde is zero outside that
! contour, so the rest of the solution can be found by solving the
! exterior Modified Helmholtz equation with just the ntilde source.
    pmaxup=0.7
    psimaxd=psimax
    wjmod=wj
    write(*,*)'psimax  psimaxd  USOL(0)    cg    rho(0)  wjmod     xk   delta '
    do k=1,1
! Calculate the interior Helmholtz coefficients
    alpha=1/2.   ! Make sure we use the correct alpha value.
    cg=-psimax**2/(2.*(psimax+wj)**2)
    xk=sqrt(-1.-2*cg)
    a2=-1.5        ! Relative amplitude of p2
    
! Express the potential on the grid via harmonics.
! Use ntilde(phi) function to deduce ntilde (rhominus).
    do j=0,nz
       zij=z(j)
       do i=0,nr
          rij=r(i)
          rsp=sqrt(zij**2+rij**2+1.e-15)
          xkrsp=xk*rsp
          if(xkrsp.lt.3*3.1415926/2)then  
             ctheta=zij/rsp
             p0=1.
             p2=(3.*ctheta**2-1.)/2
             sj0=sin(xkrsp)/xkrsp
             sj2=-sj0+(sj0-cos(xkrsp))*3./xkrsp**2
             if(xkrsp.lt.1.e-6)then;sj0=1.;sj2=0.;endif
             phirz=(psimax+wj*(1+xk**2))*(p0*sj0+a2*p2*sj2) &
                  -wj*(1+xk**2)
!             phirz=-wj+psimaxd*(p0*sj0+a2*p2*sj2)
             rhominusnew(i,j)=2.*cg*max(phirz+wjmod,0.)
          else      ! Beyond validity of internal expression. 
             rhominusnew(i,j)=0.
          endif
       enddo
    enddo
    delta=maxval(abs(rhominusnew-rhominus))
    rhominus=rhominusnew
! Solve the Modified Helmholtz equation to find the full-domain potential.    
! Charge density for phi<-Wj is provided by Helmholtz term.
    call sepeli(INTL,IORDER,rmin,rmax,nr,MBDCND,BDA,ASEPX,BDB,BSEPX,  &
          zmin,zmax,nz,NBDeli,BDC,GSEPY,BDD,XNU,                       &
          cofx,cofy,rhominus,USOL,nr+1,WORK,PRTRB,IERROR)
    if(IERROR.ne.0)write(*,*)'sepeli ERROR',IERROR,nwork,WORK(1)
    phi=0.99*USOL
    write(*,'(10f8.4)')psimax,psimaxd,USOL(0,0),cg,-rhominus(0,0),wjmod,xk,delta
    call potlntilde
    enddo
    if(.false.)then
    call autocontour(USOL,nr+1,nr+1,nz+1)
    call axlabels('r/rmax','z/zmax')
    call boxtitle('Initial Potential Contours')
    call pltend
    call autocontour(-rhominus,nr+1,nr+1,nz+1)
    call axlabels('r/rmax','z/zmax')
    call boxtitle('Initial Charge Contours')
    call pltend
 endif
 stop 
  end subroutine initntilde
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine initphi
    ar=8
    az=5
    do i=0,nr
       do j=0,nz
          USOL(i,j)=psimax*exp(-(r(i)/ar)**2-(z(j)/az)**2)
       enddo
    enddo
    phi=0.98*USOL
    psimaxd=psimax ! psimax is desired psi, psimaxd is psi in rho to get it.
  end subroutine initphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is the call that solves for USOL from rhominus:
!     call sepeli(INTL,IORDER,rmin,rmax,nr,MBDCND,BDA,ASEPX,BDB,BSEPX,  &
!          zmin,zmax,nz,NBDeli,BDC,GSEPY,BDD,XNU,                       &
!          cofx,cofy,rhominus,USOL,nr+1,WORK,PRTRB,IERROR)  
!     if(IERROR.ne.0)write(*,*)IERROR,nwork,WORK(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routines for sepeli coefficients.
  subroutine cofx(x,AFUN,BFUN,CFUN)
! Return the coefficients of the finite differences at position r=x  
! The equation form is 
!   AF(X)*D^2U/DX^2 + BF(X)*DU/DX + CF(X)*U +
!   DF(Y)*D^2U/DY^2 + EF(Y)*DU/DY + FF(Y)*U= G(X,Y)
! In a cylindrical x-domain we need (1/r)d/dr(r dp/dr)-p/dlength. So formally
! AF=1, BF=(1/r) for Laplacian.
! A Helmholtz term CF=-1./dlength might violate Diagonal Dominance.
    AFUN=1.
    BFUN=1./(x+1.e-5*rmax)     ! assumes rmin=0. avoiding division by zero.
! I suspect this is inadequate at r=0 and I think it needs examination.
    CFUN=-1./dlength**2        ! Modified Helmholtz.
  end subroutine cofx

  subroutine cofy(y,AFUN,BFUN,CFUN)
    AFUN=1.
    BFUN=0.           ! Cartesian in this direction.
    CFUN=0            ! No Modified Helmholtz term. It is in cofx.
  end subroutine cofy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine updaterhoprof(gfac)
  ! Update the negative density arising from trapped deficit
  ! in response to the phi solution updated by factor gfac. 
  ! Using the profiled wjr and psidr when lrprof=.true.
    delta=maxval(abs(USOL-phi))
    phi=(1-gfac)*phi+gfac*USOL
    if(lrprof)then  ! r-profiled trapped function.
       do i=0,nr
          wjmod=wjr(i)+psidr(i)-phi(i,0)  ! wjmod+psi=wj+psimaxd
          wjmodr(i)=wjmod
          coef=0.5*(alpha+1.5)*psidr(i)**2    &
               /max(psidr(i)+wjr(i),1.e-7)**(1.5+alpha)
          do j=0,nz
             rhominusnew(i,j)=-coef*max(phi(i,j)+wjmod,0.)**(0.5+alpha)*Fperp(i)
! Add to \tilde n the difference between Flattrap and Boltzmann 
! to use finite um, phi. Then in the Helmholtz solver
! nabla^2 phi -\phi/dl^2 = \tilde n +(n_f -\phi/dl^2 -1). 
             rhominusnew(i,j)=rhominusnew(i,j)  &
                  +(-phi(i,j)/dlength**2-1+refdenapprox(phi(i,j),um))
!             rhominusnew(i,j)=rhominusnew(i,j)  ! Without this correction.
          enddo
       enddo
       wjmod=wjmodr(0)
!     write(*,*)psidr(0),phi(0,0),wjr(0),wjmod,wjmod+phi(0,0)-wjr(0)-psidr(0)
    else  ! Unprofiled pure rho(phi).
       wjmod=wj+psimaxd-phi(0,0)  ! wjmod+psi=wj+psimaxd
       coef=0.5*(alpha+1.5)*psimaxd**2/(phi(0,0)+wjmod)**(1.5+alpha) !2CG linear
       coef=coef
       rhominusnew=-coef*max(phi+wjmod,0.)**(0.5+alpha)
!  Add on the difference between Flattrap and Boltzmann
       do i=0,nr
          do j=0,nz
             rhominusnew(i,j)=rhominusnew(i,j)  &
                  -phi(i,j)/dlength**2-1+refdenapprox(phi(i,j),um)
          enddo
       enddo
    endif
! Relaxed version of \tilde n =-2CG (phi+Wj)^(0.5+alpha))
  end subroutine updaterhoprof
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine resultwrite
    open(12,file='helmoutput.dat',status='unknown',form='unformatted')
    write(12)nz,nr
    write(12)zmax,rmax,dz,dr
    write(12)(z(i),i=0,nz),(r(i),i=0,nr)
    write(12)(fpsir(i),i=0,nr),(wjr(i),i=0,nr)
    write(12)((phi(i,j),i=0,nz),j=0,nr)
    write(12)psi0r
    close(12)
  end subroutine resultwrite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine resultread
    open(12,file='helmoutput.dat',status='old',form='unformatted')
    read(12)nzf,nrf
    if(nzf.ne.nz.or.nrf.ne.nr)then
       write(*,*)'nz,nr',nz,nr,' not correct'
       stop
    else
       read(12)zmax,rmax,dz,dy
       read(12)(z(i),i=0,nz),(r(i),i=0,nr)
       read(12)(fpsir(i),i=0,nr),(wjr(i),i=0,nr)
       read(12)((phi(i,j),i=0,nz),j=0,nr)
       read(12)psi0r
    endif
    close(12)
  end subroutine resultread
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine potlntilde
! Contours of Potential and n-tilde       
    call charsize(0.018,0.018)
    call pltinaspect(0.,rmax,0.,zmax)
    gradlx=-0.25*max(1.,zmax/rmax)
    icsw=1+16+32
    zclv(1)=1.1*rhominus(0,0);zclv(2)=0.;iclv=2 !color scaling
    call blueredgreenwhite
    call contourl(rhominus,pwork,nr+1,nr+1,nz+1,zclv,iclv,r,z,icsw)
    call gradlegend(zclv(1),zclv(2),gradlx,0.,gradlx,1.,.03,.false.)
    call axis; call axis2
    call axlabels(cyl,'z')
    call boxtitle('!Af!@ (lines)  !p!o~!o!qn (colors)')
    ilmax=nint(alog10(phi(0,0)))
    do iclv=1,nclv
       zclv(iclv)=10.**(ilmax+1)*10**(-float((iclv-1)/5))*(1.-0.2*mod(iclv+4,5))
       if(iclv.eq.nclv)exit
    enddo
    icsw=1
    call contourl(phi,pwork,nr+1,nr+1,nz+1,zclv,iclv,r,z,icsw)
    if(.false.)then ! contour rhominus
!       iclv=-1     ! Zero contour only.
!       zclv(1)=-1.e-3*rhominus(0,0)
       call color(ipurple())
       call contourl(-rhominus,pwork,nr+1,nr+1,nz+1,zclv,iclv,r,z,icsw)
       call color(15)
    endif
    call pltend
    call multiframe(0,0,0)  ! Cancel frame aspect.
  end subroutine potlntilde
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine resultplots

    call potlntilde

    if(lrprof)then
! r-profiles at z=0.       
       do j=1,nr-1
          Erbyphi(j)=0.5*(phi(j-1,0)-phi(j+1,0))/dr/phi(j,0)
       end do
       Erbyphi(0)=0.;   Erbyphi(nr)=Erbyphi(nr-1)
       call charsize(0.018,0.018)
       call pltinit(0.,rmax,0.,1.2)
       call axis; call axis2
       call axlabels(cyl,'')
       call winset(.true.)
       call color(1)
       iy=int(3*nr/4.)
       call jdrwstr(wx2nx(r(iy)),wy2ny(Erbyphi(iy))-.02,  &
            'E!d'//cyl//'!d/!af!@',1.)
       call polyline(r,Erbyphi,nr)
       call color(2)
       call polyline(r,-wjr/phi(:,0),nr+1)
       call jdrwstr(wx2nx(r(iy)),wy2ny(-wjr(iy)/phi(iy,:))+0.03, &
            '-W!d!Bj!@!d/!Ay!@('//cyl//')',1.)
       call color(3)
       call polyline(r,-fpsir,nr+1)
       iy=int(nr/10.)
       call jdrwstr(wx2nx(r(iy)),wy2ny(-fpsir(iy))-.03, &
            '|!p!o~!o!qf!d!Ay!@!d|',1.)
       call color(4)
       call polyline(r,phi(:,0)/phi(0,0),nr+1)
       call dashset(2)
!       write(*,*)psi0zy(0:3)
       call polyline(r,psi0r/psi0r(0),nr+1)
       call dashset(4)
       iy=int(nr*0.25)
       call polyline(r,phi(:,0)/psi0r,nr+1)
       call jdrwstr(wx2nx(r(iy))+.01,wy2ny(phi(iy,0)/psi0r(iy))+.02,   &
            ' !Ay!@(r)/!Ay!@!dd!d(r)',1.)
       call jdrwstr(wx2nx(r(iy)),wy2ny(phi(iy,0)/phi(0,0))-0.02,&
            '!Ay!@('//cyl//')/!Ay!@(0)',-1.05)
       call dashset(0)
       call pltend
    endif
       
       
! z-shape variation plots
    call pltinit(0.,zmax,0.,1.)
    call charsize(.018,.018)
    if(lzlog)then
       zloc=zmax*.7
       call scalewn(0,zmax,phi(0,nz)/phi(0,0),1.,.false.,.true.)
       call dashset(1)
       call polyline([0.,zloc],[1.,exp(-zloc)],2)
       call dashset(0)
       call jdrwstr(wx2nx(zloc),wy2ny(exp(-zloc)),'exp(-z)',-1.)
    endif
    call axis; call axis2
    call axlabels('z','!Af!@(z)/!Ay!@')
    call legendline(.7,.9,258,'Labels: '//cyl)
    call winset(.true.)
    nj=int(rmax/2)+1
    psi=phi(0,0)            
    do j=1,nj
       k=(j-1)*(nr/(nj-1))
       k=nint(0.4+2*(j-1)*(nr-0.5)/rmax)
       call fwrite(r(k),iwidth,0,arg)
       call color(j)
       call labeline(z,phi(k,:)/phi(k,0),nz+1,arg,iwidth)
       call dashset(j)
       if(j.eq.1)then
          zk=sqrt((psi0r(0)/(psi0r(0)+wjr(0)))**2-1.)
          hphi=1./(psi0r(0)/(psi0r(0)+wjr(0))+1.)
          zj=acos(1-(psi0r(0)+wjr(0))/(phi(0,0)*hphi))/zk
          jj=min(nz,nint(nz*zj/zmax))
          pfit=(1.-(psi0r(0)+wjr(0))/psi)*exp(-(z-zj))
          if(alpha.ne.0.5)then
               write(*,*)'WARNING: Wrong Alpha for analytic phi form',alpha
          else
             call polymark(z,1.-hphi*(1.-cos(zk*z)),jj,10)
             call polymark(z(jj),pfit(jj),nz+1-jj,11)
          endif
       endif
       call dashset(0)
    enddo
    call color(15)
    call pltend

  end subroutine resultplots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine EzErfill
    do i=1,nz-1
       do j=1,nr-1
          Ez(j,i)=-0.5*(phi(j,i+1)-phi(j,i-1))/dz
          Er(j,i)=-0.5*(phi(j+1,i)-phi(j-1,i))/dr
       enddo
    enddo
    Ez(:,0)=0.
    Ez(:,nz)=Ez(:,nz-1)
    Ez(0,:)=Ez(1,:)
    Ez(nr,:)=Ez(nr-1,:)
    Er(0,:)=-Er(1,:)
    Er(nr,:)=Er(nr-1,:)
    Er(:,0)=Er(:,1)
    Er(:,nz)=Er(:,nz-1)
  end subroutine EzErfill
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine Eplots
    izm=int(nz*0.75)
    zm=z(izm)
    call multiframe(2,1,1)
    call pltinaspect(0.,rmax,0.,zm)
    call charsize(0.018,0.018)
    call axis; call axis2
    call axlabels('','z')
    call legendline(0.5,0.9,258,'E!d!A|!@!d contours')
    call minmax2(Ez,nr,nr,izm,emin,emax)
    ilmax=nint(alog10(emax)+0.5)
    do iclv=1,nclv
       zclv(iclv)=10.**(ilmax)*10**(-float((iclv-1)/5))*(1.-0.2*mod(iclv+4,5))
       if(iclv.eq.nclv)exit
    enddo
    icsw=1
    if(lsparselabels)then
       do i=1,11    ! Reduced from iclv
          ic=-1
          if(mod(i+4,5).eq.0.or.mod(i+4,5).eq.3)ic=1
          call contourl(Ez,work,nr+1,nr+1,izm+1,zclv(i),ic,r,z,icsw)
       enddo
    else
       call contourl(Ez,work,nr+1,nr+1,izm+1,zclv,iclv,r,z,icsw)
    endif
    call pltinaspect(0.,rmax,0.,zm)
    call charsize(0.018,0.018)
    call axis; call axis2
    call axlabels(cyl,'z')
    call legendline(0.5,0.9,258,'E!d!A`!@!d contours')
    call minmax2(Er,nr,nr,izm,emin,emax)
    ilmax=nint(alog10(emax)+0.5)
    do iclv=1,nclv
       zclv(iclv)=10.**(ilmax)*10**(-float((iclv-1)/5))*(1.-0.2*mod(iclv+4,5))
       if(iclv.eq.nclv)exit
    enddo
    icsw=1
    if(lsparselabels)then
       do i=1,11    ! iclv
          ic=-1
          if(mod(i+4,5).eq.0.or.mod(i+4,5).eq.3)ic=1
          call contourl(Er,work,nr+1,nr+1,izm+1,zclv(i),ic,r,z,icsw)
       enddo
    else
       call contourl(Er,work,nr+1,nr+1,izm+1,zclv,iclv,r,z,icsw)
    endif
    call pltend
    call multiframe(0,0,0)  !Reset plot box.

    do j=1,nr
       Ezp(j)=maxval(Ez(j,:))
       Erp(j)=maxval(Er(j,:))
       Erzr(j)=Erp(j)/Ezp(j)
    enddo
    yleg=.5
    call pltinit(0.,r(nr),0.,Ezp(1))
    call charsize(0.018,0.018)
    call axis
    call polyline(r(1),Ezp(1),nr)
    call legendline(.02,yleg,0,' E!d!A|!@!d')
    call axlabels(cyl,'E!dpeak!d')
    call dashset(1)
    call polyline(r(1),Erp(1),nr)
    call legendline(.02,yleg-.05,0,' E!d!A`!@!d')
    call dashset(2)
    call scalewn(0.,r(nr),0.,1.1*maxval(Erzr(1:nr)),.false.,.false.)
    call axptset(1.,1.)
    call ticrev
    call altyaxis(1.,1.)
    call axlabels('','E!d!A`!@peak!d/E!d!A|!@peak!d')
    call ticrev
    call axptset(0.,0.)
    call polyline(r(1),Erzr(1),nr)
    call legendline(.7,yleg,0,' E!d!A`!@!d/E!d!A|!@!d')    
    call dashset(0)
    call pltend
  end subroutine Eplots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine outputphi(nz,ny,z,y,phiguarded,phiyz)
    integer nz,ny
    real phiguarded(-1:nz+1,-1:ny+1)
    real phiyz(0:ny,0:nz)
    real z(0:nz),y(0:ny)
! phi(0:nr,0:nz) plays the prior role of phizy so no longer
! phiguarded(0:nz,0:ny)=phizy, instead transpose.
    do i=0,ny
       do j=0,nz
          phiguarded(j,i)=phiyz(i,j)
       enddo
    enddo
! Fill guard values
    phiguarded(-1,:)=phiguarded(1,:) ! z0=0 symmetry
    phiguarded(:,-1)=phiguarded(:,2) ! y0=-dy/2 symmetry
    phiguarded(nz+1,:)=2.*phiguarded(nz,:)-phiguarded(nz-1,:) 
    phiguarded(:,ny+1)=2.*phiguarded(:,ny)-phiguarded(:,ny-1)
! Output.
    open(14,file='helmphiguard.dat',status='unknown',form='unformatted')
    write(14)nz,ny
    write(14)z,y
    write(14)phiguarded
!    write(14)((phiguarded(i,j),i=-1,nz+1),j=-1,ny+1)
    close(14)
  end subroutine outputphi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine inputphi(Lz,nzf,nyf,z,y,phiguarded)
! This reads into phiguarded as if it were 
!                 real phiguarded(-1:nzf+1,-1:nyf+1)
! Subsequently access phiguarded accordingly. It is like phi transposed.
! This routine is not used here, but it is a template. If output is
! changed, the version in (e.g.) orbitpoincare needs to be updated.
    integer Lz,nzf,nyf
    real phiguarded(*)
    real z(0:Lz),y(0:Lz)
    open(15,file='helmphiguard.dat',status='old',form='unformatted')
    read(15)nzf,nyf
    if(nzf.gt.Lz.or.(nzf+3)*(nyf+3).gt.(Lz+3)**2) &
        stop 'Incompatible phiguarded ranks'
    write(*,*)nzf,nyf
    read(15)(z(i),i=0,nzf),(y(i),i=0,nyf)
    read(15)((phiguarded(i+(j-1)*(nzf+3)),i=1,nzf+3),j=1,nyf+3)
    close(15)
    write(*,*)'Read in ftfphi.dat successfully'
  end subroutine inputphi
!*********************************************************************
end module helmhole
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*********************************************************************
! Pade approximation to Z-function for flattrap density.
real function refdenapprox(phi,um)
  real, parameter :: sqrt_pi = 1.77245385
  d1=.829
  d2=.363
  slope=(1.-(2-d1)*um**2)/(1+d1*um**2+d2*um**4)
  umfac=(1.+1.2*um**2)/((1+.2*um**2+2.*um**6)*2./sqrt_pi)
  sqphi=sqrt(max(phi,0.))
!      refdenapprox=1.+  slope*phi                                        &
!     &     /(1+.895*slope*umfac*sqphi)
! Correction for substantial um, a modification of the above comment
  refdenapprox=1.+ (slope+0.13*min(1.5,um**8)*sqphi)*phi         &
       &     /(1+.895*slope*umfac*sqphi)
end function refdenapprox
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function Vflatapprox(phi,um)
! Approximation to the Classical Potential for a flattrapped distribution
! in a moving Maxwellian, using the above rational approximation.
! (Without the substantial um correction)
  real, parameter :: sqrt_pi = 1.77245385
  d1=.829
  d2=.363
  slope=(1.-(2-d1)*um**2)/(1+d1*um**2+d2*um**4)
  umfac=(1.+1.2*um**2)/((1+.2*um**2+2.*um**6)*2./sqrt_pi)
  xfac=0.895*slope*umfac
! Now we must integrate x/(1+sqrt(x)) which is 
! (2x^3/2/3-x+2x^1/2-2\ln(1.+x^1/2))
! But that only works for positive xfac
!      if(slope*phi.lt.4.e-3.or.xfac.lt.0.)then
  if(xfac.lt.0.)then
     Vflatapprox=slope*phi**2/2.
  else
     sx=xfac*sqrt(phi)
     if(sx.lt.0.1)then
! Avoid rounding errors giving inaccurate cancellation in the full expression
! Expand log(1+sx) = sx-sx^2/2+sx^3/3-sx^4/4+sx^5/5 ...
        Vflatapprox=2.*(sx**4/4.-sx**5/5.+sx**6/6.-sx**7/7.)*slope/xfac**4
     else
        Vflatapprox=(sx*(2.-sx*(1-2.*sx/3.))-2.*alog(1.+sx))*slope/xfac**4
     endif
  endif
  Vflatapprox=-Vflatapprox
end function Vflatapprox
!*********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zprofiles
! Run with various parameter settings alpha and Wj and plot z-profiles.
  use helmhole
  character*100 thelabel(nprof)
  mprof=3
  do kc=1,mprof
     if(kc.eq.1)then; psimax=0.1;alpha=.02;Wj=-.028; endif
     if(kc.eq.2)then; psimax=0.1;alpha=1;Wj=-.02; endif
     if(kc.eq.3)then; psimax=0.1;alpha=.5;Wj=-.002; endif
     write(thelabel(kc),101)alpha,-Wj
101  format(' !Aa!@=',f4.2,', -W!Bj!@=',f5.3)
     call initrzbd
     call initwbpsi
     gfacp=1.
     write(*,*)' rhominus(0)      phi(0)     psimaxd        wjmod       delta'
     do k=1,niter
        call updaterhoprof(gfacp)
        psimaxd=psimaxd+psiupfac*(psimax-USOL(0,0)) ! Update dynamic psimax
        psidr=psidr+psiupfac*(psi0r-USOL(:,0))
        rhominus=rhominusnew
        if(k.le.50.or.k.eq.niter)write(*,'(i4,6e13.5)') &
             k,rhominus(0,0),USOL(0,0) &
             ,psimaxd,wjmod,delta
        call sepeli(INTL,IORDER,rmin,rmax,nr,MBDCND,BDA,ASEPX,BDB,BSEPX,  &
             zmin,zmax,nz,NBDeli,BDC,GSEPY,BDD,XNU,                       &
             cofx,cofy,rhominus,USOL,nr+1,WORK,PRTRB,IERROR)  
        if(IERROR.ne.0)write(*,*)'sepeli ERROR',IERROR,nwork,WORK(1)
        if(delta.lt.1.e-6*psi0r(0))exit
     enddo
     phinprof(:,kc)=phi(0,:)
     rhonprof(:,kc)=-rhominus(0,:)
  enddo
  call pfset(ips)
  call multiframe(2,1,1)
  call pltinit(0,zmax,0.,1.1*psimax)
  call charsize(0.016,0.016)
  call axis; call axis2
  call axlabels('','   !Af!@(z)')
  do kc=1,mprof
     call color(kc)
     call dashset(kc-1)
     call polyline(z,phinprof(0,kc),nz+1)
     call legendline(0.5,1.-0.08*kc,0,thelabel(kc))
  enddo
  call color(15)
  call pltinit(0,zmax,0.,rhonprof(0,2))
  call charsize(0.016,0.016)
  call axis; call axis2
  call axlabels('z','|!p!o~!o!qn(z)|    ')
  do kc=1,mprof
     call color(kc)
     call dashset(kc-1)
     call polyline(z,rhonprof(0,kc),nz+1)
     call legendline(0.5,1.-0.08*kc,0,thelabel(kc))
  enddo
  call pltend
  
end subroutine zprofiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mainhelmhole
! Iterates a solution to a problem where the extra source of charge
! on the RHS of Helmholtz equation is a function of potential only.
! It demonstrates that (1) The charge function must be referenced to
! a fixed peak target potential (updaterho2) not to zero (updaterho)
! (2) That initial anisotropic potential profile relaxes to cylindrical
! symmetry when rho has no built-in radial (anisotropic) dependence.
  use helmhole

  call parsearguments
  if(lzprof)then  ! Special plot of different z-profiles hard coded
     call zprofiles
  else
  call initrzbd   ! Default case, document equilibrium set by command line.
  if(lrprof)then  ! Use shaped transverse profiles
     call initwbpsi
  else           ! Use global n(phi) and initial gaussians.
     write(*,*)'Calling initntilde'
     call initntilde
!     call initphi
     psi0r=psimax
     wjr=wj
  endif
  
  call pltinit(0,rmax,-.01*psimax,(1.+alpha+abs(wj)/psimax)*psimax)
  call charsize(0.018,0.018)
  call axis
  call boxtitle('Iterations to convergence')
  call axlabels('r','phi (line)     rho (dashed) at z=0')
  gfacp=1.+0.5/dlength**2
  gfacp=1.  ! This works better for anisotropic, no dynamic
  write(*,*)' rhominus(0)      phi(0)     psimaxd        wjmod       delta'
  do k=1,niter
     call updaterhoprof(gfacp)
     psimaxd=psimaxd+psiupfac*(psimax-USOL(0,0)) ! Update dynamic psimax
     psidr=psidr+psiupfac*(psi0r-USOL(:,0))
     rhominus=rhominusnew
     
     if(k.le.50.or.k.eq.niter)write(*,'(i4,6e13.5)')k,rhominus(0,0),USOL(0,0) &
          ,psimaxd,wjmod,delta
     call color(mod(k,14)+1)
     call polyline(r,USOL,nr+1)
     call dashset(2)
     call polyline(r,-rhominus,nr+1)
     call dashset(0)

     call sepeli(INTL,IORDER,rmin,rmax,nr,MBDCND,BDA,ASEPX,BDB,BSEPX,  &
          zmin,zmax,nz,NBDeli,BDC,GSEPY,BDD,XNU,                       &
          cofx,cofy,rhominus,USOL,nr+1,WORK,PRTRB,IERROR)  
     if(IERROR.ne.0)write(*,*)'sepeli ERROR',IERROR,nwork,WORK(1)
     if(delta.lt.1.e-6*psi0r(0))exit
  enddo
  call pltend

  call resultplots

  call EzErfill
  call Eplots
  call outputphi(nz,nr,z,r,phiguarded,phi)

endif
  
end subroutine mainhelmhole
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program hole
  call mainhelmhole
end program hole
