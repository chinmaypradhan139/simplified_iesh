module mod_global

implicit none
  integer::nsize,INFO,total_dimensions,Hi,Ne,nold!(be careful with setting dt) 
  real*8 :: tau,planks_constant,kT,Vr,Energy_of_diabat,U0,dt,dtc,total_time
  complex*16 :: population
  real*8 :: time,rnd,omega,gh,dG,Band_width,RAND,rhos
  complex*16,dimension(:,:), allocatable :: c
  real*8,dimension(:,:), allocatable :: Energy_hamil,g,Identity,H,Grd,vdij
  real*8,dimension(:,:,:), allocatable :: Gradient_hamil,acw
  real*8,dimension(:), allocatable :: Energy,acc1,E_levels,w_levels,knot_x
  real*8,dimension(:), allocatable ::pos,v,old_acc2,momentum,population_mat
  real*8, dimension(:), allocatable :: mass,pop_mat
  integer,dimension(:), allocatable :: lambda
  complex*16, dimension(:), allocatable :: cwork
  integer(kind=4) ::iseed
  integer :: ham_sp
  integer :: n_q
  

contains
!.......................................................................
subroutine setup_initial_values1
implicit none
integer(kind=4) :: O,j
integer(kind=4), dimension(:),allocatable::x
real*8 :: wt
integer:: ut,xt,ip,kn
real*8, dimension(14) :: inpot

 open(25,file='fort.23')
 do ip=1,14
 read(25,*) inpot(ip)
 enddo
 close(25)

write(11,*) inpot

Hi=inpot(1)
Ne=int(Hi/2)
tau=inpot(4)
Band_width=inpot(5)
!mass=inpot(6)
omega=inpot(7)
gh=inpot(8)
dG=inpot(9)
KT=inpot(10)
dtc=inpot(11)
dt=inpot(12)
total_dimensions=inpot(13)
wt=inpot(14)
total_time=wt*50
rhos=real(Hi)/Band_width
Vr=sqrt(tau/(2*3.1416))
ham_sp=1
  allocate(pos(total_dimensions))
  allocate(v(total_dimensions))
  allocate(momentum(total_dimensions))
  allocate(acc1(total_dimensions))
  allocate(old_acc2(total_dimensions))
  allocate(acw(Hi,Hi,total_dimensions))
  allocate(H(Hi,Hi))
  allocate(Energy_hamil(Hi,Hi)) 
  allocate(Energy(Hi))
  allocate(Gradient_hamil(Hi,Hi,total_dimensions))
  allocate(c(Ne,Hi))
!  allocate(b(Ne,Hi))
!  allocate(A(Ne,Hi))
  allocate(g(Ne,Hi))
  allocate(lambda(Ne))
  allocate(Identity(Hi,Hi))
  allocate(Grd(Hi,Hi))
  allocate(vdij(Hi,Hi))
  allocate(population_mat(int(total_time/dtc)))
  allocate(E_levels(int(Hi/2)))
  allocate(knot_x(int(Hi/2)))
  allocate(w_levels(int(Hi/2)))
  allocate(mass(total_dimensions))
  allocate(pop_mat(int(total_time)))
call random_seed(size=O)
  allocate(x(O))
   do j=1,O
   x(j)=j**6*iseed+2777772
   enddo
  call random_seed(put=x)
 open(2,file='rndom',status='new')
 write(2,*) iseed,x,O
 


open(167,file='raw_x.txt')
do kn=1,int(Hi/2)
read(167,*) knot_x(kn)
enddo


open(169,file='raw_w.txt')
do xt=1,int(Hi/2)
read(169,*)w_levels(xt)
enddo
close(169)





mass(1)=inpot(6)

  

end subroutine
!.........................................................
  
subroutine gaussian_random_number(rnd0)
!USE IFPORT
   !! generates gaussian distribution with center 0, sigma 1
   !! q0+sig*rnd gives center=q0, sigma=sig
   implicit none
   integer(kind=4) :: n,j,M,O,k
   real*8,intent(out)::rnd0
   real*8 rnd1,rnd2,pi
   pi=dacos(-1.d0)
   call random_number(rnd1)
   call random_number(rnd2)
   rnd0 = dsqrt(-2*log(rnd1))*dcos(2*pi*rnd2)

end subroutine gaussian_random_number
!.............................................................
subroutine setup_initial_values2
integer :: n,m,p,y,i
real*8 :: rnd2,rnd1

call gaussian_random_number(rnd2)
call gaussian_random_number(rnd1)
momentum(1)=sqrt(mass(1)*KT)*rnd2
v(1)=(momentum(1)/mass(1))  

pos(1)=rnd1/(sqrt(mass(1)*omega**2/KT))
time=0.00




c=0
do i=1,Ne
   c(i,i+1)=1
enddo 














do n=1,Ne
lambda(n)=n
enddo





call potential(ham_sp)

c=matmul(c,Energy_hamil)



end subroutine 
!........................................................................
subroutine diag_wrapper(matrix,nsize,eigen_values,eigen_vectors)
real*8, intent(inout) :: matrix(nsize,nsize)
real*8, dimension(nsize,nsize) :: mat
integer LWORK,nsize
real*8, allocatable :: WORK(:)
real*8, intent(out) :: eigen_vectors(nsize,nsize),eigen_values(nsize)
mat=matrix
LWORK=3*nsize-1
allocate(WORK(LWORK))
call dsyev('V','U',nsize,mat,nsize,eigen_values,WORK,LWORK,INFO)
eigen_vectors=mat

end subroutine

!..........................................................................
subroutine logm(mat,log_mat,n)
   !! http://arxiv.org/pdf/1203.6151v4.pdf
   implicit none
   integer,intent(in):: n
   real*8,intent(in):: mat(n,n)
   real*8,intent(out):: log_mat(n,n)
   integer i
   complex*16 T(n,n),en(n),vect(n,n)
   complex*16 dd(n,n)
    
   call schur(mat,T,n,en,vect,nold,cwork)
   dd=0.d0
   do i=1,n
     dd(i,i)=cdlog(t(i,i)/cdabs(t(i,i)))
   enddo
      
   log_mat=matmul(vect,matmul(dd,conjg(transpose(vect))))

end subroutine logm


!.............................................................................

subroutine schur(mat,T,n,eigen_value,eigen_vect,nold,cwork)
   !! Diaganalizing matrix using dsyevr. First m_values eigen values and
!eigenvectors computed.
   !! The module's common variables should contain:

   !! Initialize nold=0

   !! nold makes sure that everytime value of n changes, work and iwork
!are re-allocated for optimal performance.
   !! mat is destroyed after use.

   implicit none
   integer,intent(in) :: n
   integer,intent(inout) :: nold
   complex*16,intent(out) :: eigen_value(n),eigen_vect(n,n)
   real*8,intent(in) :: mat(n,n)
   complex*16,intent(out) :: T(n,n)
   complex*16,allocatable,intent(inout):: cwork(:)
   real*8 rwork(n)
   complex*16 mat_c(n,n)

   integer lwork
   logical:: select
   logical bwork(n)
   integer sdim,info,AllocateStatus

   T=mat

   info=0
   sdim=0

   if(nold.ne.n .or. .not.allocated(cwork)) then
   !if(nold.ne.n) then
     lwork=-1
     if(allocated(cwork))deallocate(cwork)
     allocate(cwork(n))
     call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)
     lwork=int(cwork(1))
     deallocate(cwork)
     allocate(cwork(lwork),STAT=AllocateStatus)
     if(allocatestatus.ne.0) write(6,*)"problem in schur, allocation"
     nold=n
   endif

   lwork=size(cwork)
   call zgees('V','N',SELECT,N,T,n,SDIM,eigen_value,eigen_vect,n,cWORK,LWORK,rwork,BWORK,INFO)

 !  write(57,*) mat
   if(info.ne.0) then
     write(6,*) "problem in scur",info
     stop
   endif

end subroutine schur
!...................................................................................................................
subroutine potential(sp)
implicit none
real*8 :: U1
integer :: i,j,k,l,x,y,r
integer, intent(in) :: sp

If (sp.eq.1) then    
    U0=0.5*mass(1)*(omega**2)*(pos(1))**2
    U1=0.5*mass(1)*(omega**2)*(pos(1)-gh)**2+dG
    H(1,1)=U1-U0
    do i=2,Hi
        if (i.le.int(Hi/2)) then
                H(1,i)=sqrt(Band_width*w_levels(i))*Vr/2
                H(i,1)=sqrt(Band_width*w_levels(i))*Vr/2
        else
                H(1,i)=sqrt(Band_width*w_levels(i-int(Hi/2)))*Vr/2
                H(i,1)=sqrt(Band_width*w_levels(i-int(Hi/2)))*Vr/2
        end if
        do j=2,Hi
                if (i.eq.j) then
                        if (i.le.(int(Hi/2)+1)) then
                                H(i,j)=-(Band_width/2)*(0.5+0.5*knot_x(int(Hi/2)-i+2))
                        else
                                H(i,j)=(Band_width/2)*(0.5+0.5*knot_x(i-int(Hi/2)-1))
                        end if
                else
                        H(i,j)=0.0
                end if
        end do
   enddo

else
    Vr=sqrt(tau/(2*3.1416*rhos))
    U0=0.5*mass(1)*(omega**2)*(pos(1))**2
    U1=0.5*mass(1)*(omega**2)*(pos(1)-gh)**2+dG
    H(1,1)=U1-U0
    do j=2,Hi
      do i=2,Hi
        if (i.eq.j) then
          H(i,j)=i*(Band_width/(Hi-2))-Band_width*(Hi+2)/(2*Hi-4)
        else
          H(i,j)=0
        endif
      enddo
    enddo
    do r=2,Hi
    H(r,1)=Vr
    H(1,r)=Vr
    enddo
end if







do x=1,Hi
 do y=1,Hi
   write(217,'(f12.8$)')H(x,y)
 enddo
 write(217,*)
enddo


!if (time.ge. 593690) then
!    write(56,*) time,pos
!end if



nsize=Hi
call diag_wrapper(H,nsize,Energy,Energy_hamil)






Gradient_hamil(1,1,1)=-mass(1)*(omega**2)*gh
do k=2,Hi
do l=2,Hi
Gradient_hamil(l,1,1)=0
Gradient_hamil(1,k,1)=0
if (l.eq.k) then
Gradient_hamil(l,k,1)=0
else
Gradient_hamil(l,k,1)=0
endif
enddo
enddo





end subroutine
!...........................................................................
subroutine nonadiabaticvector(t,u)
implicit none
integer :: acwi
integer, intent(in) :: t,u

!write(94,*) lambda(t)

acw=0

do acwi=1,total_dimensions

if (t.ne.u) then
acw(lambda(t),u,acwi)=(sum(Energy_hamil(:,lambda(t))*matmul(Gradient_hamil(:,:,acwi),Energy_hamil(:,u))))/(Energy(lambda(t))-Energy(u))
else
acw(lambda(t),u,acwi)=0
endif

end do
end subroutine
!..........................................................................
subroutine force
implicit none
integer :: acc1i,z,u,q
real*8, dimension(total_dimensions,Ne) :: m_acc1
!call potential
do acc1i=1,total_dimensions

do z=1,Ne
m_acc1(acc1i,z)=-((sum((Energy_hamil(:,lambda(z)))*matmul(Gradient_hamil(:,:,acc1i),Energy_hamil(:,lambda(z)))))/mass(1))
end do
end do
acc1=sum(m_acc1(:,1:Ne))-omega**2*pos(1)
end subroutine
!..........................................................................
subroutine Rungekutta
implicit none
integer :: i,j,p,q
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M

do i=1,Hi
do j=1,Hi 
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do

do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-sum(v*acw(p,q,:))
enddo
enddo

k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))
c=c+(k1+2*k2+2*k3+k4)/6
end subroutine
!.........................................................................

subroutine Rungefutta
implicit none
integer :: i,j,p,q,x,y
complex*16,dimension(Hi,Hi) :: Vad
complex*16,dimension(Ne,Hi) ::k1,k2,k3,k4
complex*16,dimension(Hi,Hi) ::M
!real*8,intent(in), dimension(Hi,Hi) :: dij
do i=1,Hi
do j=1,Hi
if (i.eq.j) then
Vad(i,j)=Energy(i)/(cmplx(0,1))
else
Vad(i,j)=0
end if
end do
end do


do p=1,Hi
do q=1,Hi
M(p,q)=Vad(p,q)-vdij(p,q)


enddo
enddo





k1=dt*matmul(c,transpose(M))
k2=dt*matmul((c+k1/2),transpose(M))
k3=dt*matmul((c+k2/2),transpose(M))
k4=dt*matmul((c+k3),transpose(M))

c=c+(k1+2*k2+2*k3+k4)/6
end subroutine



!...................................................................................
subroutine evolve_quantum
implicit none
integer, dimension(Ne) ::p
integer :: i,d,r,s
complex*16 :: det_S1,det_S2
complex*16, dimension(Ne,Ne) ::  S1,S2
complex*16 :: bdi,Adi
real :: Aii,gid,rnd


call rungefutta

call random_number(rnd)
jloop : do i=1,Ne
  do d=1,Hi
    if (FINDLOC(lambda,d,1).eq.0) then
        p=lambda

       do r=1,Ne
         S1(:,r)=c(:,p(r))
       enddo

       call modulus(S1,Ne,det_S1)

       p(i)=d


       do s=1,Ne
         S2(:,s)=c(:,p(s))
       enddo

       call modulus(S2,Ne,det_S2)



       Adi=det_S2*conjg(det_S1)
       bdi=-2*real(Adi*vdij(d,i))
       Aii=det_S1*conjg(det_S1)
       gid=dt*real(bdi/Aii)
    
     

       if (gid>rnd) then
          write(22,*) time,'tim'
          call nonadiabaticvector(i,d)
          write(22,*) lambda
          write(22,*) lambda(i),d
!write(22,*)sum(v*acw(lambda(i),d,:))/norm2(acw(lambda(i),d,:))**2+(2*Energy(lambda(i))/mass(1))-(2*Energy(d)/mass(1))
          call hop1(i,d)
          write(22,*) lambda
          write(22,*) gid,rnd
!write(22,*) sum(v*acw(lambda(i),d,:))/norm2(acw(lambda(i),d,:))**2+(2*Energy(lambda(i))/mass(1))-(2*Energy(d)/mass(1))             

     !     write(17,*) time+real(11-n_q)*dt,lambda
          exit jloop
          write(17,*) 'a'
       endif
     end if
   enddo
enddo jloop

write(16,*) time,lambda


end subroutine
!.........................................................................
subroutine classical_evolution
integer :: p,r,TT,yt,i
real*8, dimension(Hi,Hi) :: old_Energy_hamil
real*8,dimension(Hi) :: signature
real*8 :: KE,TE,EE,rnd1
call potential(ham_sp)
call force



TT=int(total_time/dtc)

do while(time.le.total_time)


write(36,*) time,pos(1)

old_Energy_hamil=Energy_hamil
call velocity_verlet


call signt(old_Energy_hamil,signature)

do r=1,Hi
 if (signature(r)<0) then
    Energy_hamil(:,r)=-Energy_hamil(:,r)
 endif
enddo

call vdotd(old_Energy_hamil)

n_q=int(dtc/dt)
do while(n_q>0) 

call evolve_quantum

n_q=n_q-1
enddo

write(33,*) time,pos(1),v(1)

call populations
yt=int(time/dtc)
population_mat(yt)=real(population)



time=time+dtc
enddo




end subroutine
!................................................................................
subroutine hop1(t,u)
implicit none
integer, intent(in) :: t,u
real*8 :: para_v
real*8, dimension(:),allocatable :: perp_v
real*8 :: cond

write(22,*) lambda(t),u,'dd'

allocate(perp_v(total_dimensions))

cond=((sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:)))**2+(2*Energy(lambda(t))/mass(1))-(2*Energy(u)/mass(1)))

write(22,*) cond

if (cond>0) then
    write(22,*) 'cc'
    para_v=sum(v*acw(lambda(t),u,:))/norm2(acw(lambda(t),u,:))
    perp_v=v-para_v*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))
    para_v=((para_v)/abs(para_v))*sqrt(para_v**2+(2*Energy(lambda(t))/mass(1))-(2*Energy(u)/mass(1)))
    v=(para_v)*acw(lambda(t),u,:)/norm2(acw(lambda(t),u,:))+perp_v
    lambda(t)=u
end if

!Note: (sum(v*acw(lambda(t),u,:))/norm(acw(lambda(t),u,:))) is the parallel component of velocity in the direction of nonadiabatic coupling vector
end subroutine
!................................................................................


subroutine hop(t,u)
implicit none
integer, intent(in) :: t,u
integer :: reciprocal_mass_loop,v_loop
real*8, dimension(total_dimensions) :: reciprocal_mass
real*8 :: a,b,gama,frustrated_condition




do reciprocal_mass_loop=1,total_dimensions
     reciprocal_mass(reciprocal_mass_loop)=1/mass(reciprocal_mass_loop)
end do



       a=0.5*sum(reciprocal_mass*acw(lambda(t),u,:))
       b=sum(v*acw(lambda(t),u,:))
     
       write(74,*)time,acw(lambda(t),u,:),lambda(t),t,u           
   
       frustrated_condition=b**2+4*a*(Energy(lambda(t))-Energy(u))
       if (time.ge.593690) then
          write(75,*) time,frustrated_condition
       end if

       if (frustrated_condition>0) then
           if (b<0) then
           gama=(b+sqrt(frustrated_condition))/2*a
           else
           gama=(b-sqrt(frustrated_condition))/2*a
           end if
           do v_loop=1,total_dimensions
             v(v_loop)=v(v_loop)-gama*acw(lambda(t),u,v_loop)/mass(v_loop)
           end do
           lambda(t)=u
        end if




end subroutine hop
!.................................................................................................
subroutine velocity_verlet
real*8 :: delr(total_dimensions),delv(total_dimensions)
real*8 :: gama_dt,gamma_B,c0,c1,c2
gamma_B=2*omega
gama_dt=gamma_B*dtc
c0=dexp(-gama_dt)
c1=1.d0/gama_dt*(1.d0-c0)
c2=1.d0/gama_dt*(1.d0-c1)
call stochastic_force(delr,delv)
pos=pos+c1*dtc*v+c2*dtc*dtc*acc1+delr
old_acc2=acc1

!if (time.ge.593690) then
!    write(61,*) time,pos
!    write(62,*) time,v
!    write(63,*) time,c1
!    write(64,*) time,c2
!    write(65,*) time,acc1
!    write(66,*) time,delr
!end if


call potential(ham_sp)
call force
v=c0*v+(c1-c2)*dtc*old_acc2+c2*dtc*acc1+delv

!if (time.ge.593690) then
!    write(67,*) time,pos
!    write(68,*) time,v
!    write(69,*) time,c1
!    write(70,*) time,c2
!    write(71,*) time,acc1
!    write(72,*) time,delv
!    write(73,*) time,old_acc2
!end if




end subroutine
!................................................................................
subroutine vdotd(old_Energy_hamil)
real*8, intent(in),dimension(Hi,Hi) :: old_Energy_hamil
!variable for checking sign of the wave
!real*8, intent(out),dimension(Hi,Hi) :: dij
real*8, dimension(Hi,Hi) :: Ut,logarithm_Ut
integer :: i,x,y

!if (time.ge.593690) then
! do x=1,Hi
!  do y=1,Hi
!   write(217,*)time,Energy_hamil(x,y)
!  enddo
!   write(217,*)
! enddo
!end if


!if (time.ge.593690) then
! do x=1,Hi
!  do y=1,Hi
!   write(217,*)time,old_Energy_hamil(x,y)
!  enddo
!   write(217,*)
! enddo
!end if






Ut=matmul(transpose(old_Energy_hamil),Energy_hamil)
!if (time.ge.593690) then
!    write(58,*) time,Ut
!end if

call logm(Ut,logarithm_Ut,Hi)
vdij=logarithm_Ut/dtc




end subroutine
!........................................................................................
subroutine signt(old_Energy_hamil,signature)

real*8, intent(in),dimension(Hi,Hi) :: old_Energy_hamil
!variable for checking sign of the wave

integer :: i
real*8, intent(out),dimension(Hi) :: signature
do i=1,Hi
signature(i)=sum(old_Energy_hamil(:,i)*Energy_hamil(:,i))
enddo



end subroutine

!....................................................................................
subroutine velocity_verlet2
pos=pos+v*dtc+0.5*acc1*dtc*dtc
old_acc2=acc1
call potential(ham_sp)
call force
v=v+0.5*(acc1+old_acc2)*dtc
end subroutine

!...................................................................................
subroutine stochastic_force(delr,delv)
real*8, intent(out) :: delr(total_dimensions),delv(total_dimensions)
integer :: i
real*8 :: rnd1,rnd2,sig_r,sig_v,sig_rv,gdt,gamma_B
gamma_B=2*omega
gdt=gamma_B*dtc

do i=1,total_dimensions

sig_r=dtc*dsqrt(KT/mass(i)*1.d0/gdt*(2-1.d0/gdt*(3-4*dexp(-gdt)+dexp(-2*gdt))))
sig_v=dsqrt(KT/mass(i)*(1-dexp(-2*gdt)))
sig_rv=(dtc*KT/mass(i)*1.d0/gdt*(1-dexp(-gdt))**2)/(sig_r*sig_v) !correlation coefficient

call gaussian_random_number(rnd1)
call gaussian_random_number(rnd2)
delr(i)=sig_r*rnd1
delv(i)=sig_v*(sig_rv*rnd1+dsqrt(1-sig_rv**2)*rnd2)
enddo

end subroutine stochastic_force
!....................................................................................
subroutine populations
integer :: a,i,j,k
real*8, dimension(Hi,Hi,Ne) :: rho_a
real*8, dimension(Hi,Hi) :: rho_d
real*8, dimension(Ne) :: LUMO_population
do a=1,Ne
do i=1,Hi
do j=1,Hi
rho_a(i,j,a)=c(a,i)*conjg(c(a,j))
enddo
enddo
enddo

do k=1,Ne
rho_d=matmul(Energy_hamil,(matmul(rho_a(:,:,k),transpose(Energy_hamil))))
LUMO_population(k)=rho_d(1,1)
enddo


population=sum(LUMO_population(1:Ne))
write(15,*)time,1-real(population)!,lambda
pop_mat(int(time/dtc))=1-real(population)

!write(18,*) time,pop_mat(int(time/dtc))
end subroutine
!...................................................................................................................

subroutine modulus(matrix,n,determinant)
 IMPLICIT NONE
     complex*16, DIMENSION(n,n) :: matrix
     INTEGER, INTENT(IN) :: n
     complex*16 :: m, temp
     INTEGER :: i, j, k, l
     LOGICAL :: DetExists = .TRUE.
     complex*16,intent(out) :: determinant
     l = 1
     !Convert to upper triangular form
     
     DO k = 1, n-1
         IF (matrix(k,k) == 0) THEN
             DetExists = .FALSE.
             DO i = k+1, n
                 IF (matrix(i,k) /= 0) THEN
                     DO j = 1, n
                         temp = matrix(i,j)
                         matrix(i,j)= matrix(k,j)
                         matrix(k,j) = temp
                     END DO
                     DetExists = .TRUE.
                     l=-l
                     EXIT
                 ENDIF
             END DO
             IF (DetExists .EQV. .FALSE.) THEN
                 determinant= 0
                 return
             END IF
         ENDIF
         DO j = k+1, n
             m = matrix(j,k)/matrix(k,k)
             DO i = k+1, n
                 matrix(j,i) = matrix(j,i) - m*matrix(k,i)
             END DO
         END DO
     END DO

     !Calculate determinant by finding product of diagonal elements
     determinant= l
     DO i = 1, n
         determinant= determinant* matrix(i,i)
     END DO

END subroutine modulus

!......................................................................................
end module














