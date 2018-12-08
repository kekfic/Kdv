program kdvadim
integer :: i,j, N,Npassi,Ncicli,Nsalva,Jciclo,Jtempo
real, allocatable, dimension(:) :: U,F,Unext,x,Err,Xnew,vel,xx_max,vv_max,aa_max,uu_max
real:: deltax,deltat,pi,errore,x_max,Umax,tempo,err_rel
character*30 nomefile
character car0,car1,car2
integer j0,j1,j2,i_max,kmax


 pi=4.0*atan(1.0)
 
open(16,file='massimo.txt')
!---------------------------------file massimo
open(10,file='Parametri.txt')
read(10,*)N
  write(*,*)'N=',N
read(10,*)xmin
  write(*,*)'xmin=',xmin
read(10,*)xmax
  write(*,*)'xmax=',xmax
read(10,*)
  write(*,*)'-----'
read(10,*)deltat
  write(*,*)'dt=',deltat
read(10,*)Ncicli
  write(*,*)'Ncicli=',Ncicli
read(10,*)Nsalva
  write(*,*)'Nsalva=',Nsalva

close(10)


write(*,*)
write(*,*)'Lettura = .....................................[OK]'
allocate(U(0:N-1),F(0:N-1),Unext(0:N-1),x(0:N-1),Err(0:N-1),vel(1:Ncicli),Xnew(1:Ncicli))
write(*,*)'Allocazione variabili = .......................[OK]'

deltax = (xmax-xmin)/N
write(*,*)'deltax=',deltax

Npassi=Ncicli*Nsalva
write(*,*)Ncicli,Npassi,Nsalva

do j=0, N-1
x(j) = xmin + j*deltax
U(j)=2*cosh(x(j))**(-2)
end do

open(3,file='Errore.txt')!*************************aprofile per l'errore

write(*,*)'Inizializzazione = : ..........................[OK]'


!----------------------------------------------
do jciclo=1,Ncicli
!----------------------------------------------


write(*,*)jciclo,U(N/2)
!----------------------------------------------
do jtempo = 0,Nsalva-1
!----------------------------------------------

call calcolaeffe(U,F,N,deltax)

do i=0,N-1
Unext(i)=U(i)+0.5*deltat*F(i)
end do

call calcolaeffe(Unext,F,N,deltax)

do i=0, N-1
   U(i)=U(i)+deltat*F(i)    
 
end do

Umax=U(0)
do i=1,N-1
 
 if (U(i).gt.Umax) then
           i_max=i
           Umax=U(i)  
    end if
    
    
     if (jciclo==25) then     
        Err(i)=abs((U(i)-2*cosh(x(i))**(-2)))    
   end if
   
   
 end do
 
  
 
!----------------------------------------------
end do
!----------------------------------------------
!Umax=maxval(U)
  x_max = x(i_max)+0.25*deltax*(U(i_max+1)-U(i_max-1))/(U(i_max)-0.5*(U(i_max+1)+U(i_max-1)))
  !x_max=deltax*i_max
  tempo=deltat*(jtempo+(jciclo-1)*Nsalva)
   write(16,*) tempo, x_max, Umax
   write(*,*) 'imax=...................',i_max





!...definisco il nome del file (kdv001.txt kdv002.txt ecc.)
 j2=Jciclo/100
 j1=(Jciclo-100*j2)/10
 j0= Jciclo-100*j2-10*j1

 car2=char(j2+ichar('0'))
 car1=char(j1+ichar('0'))
 car0=char(j0+ichar('0'))

 
nomefile='kdv'//car2//car1//car0//'.txt'

write(*,*)'Fine ciclo numero ',jciclo,'/',Ncicli
 
write(*,*)'Inizio scrittura su file ',nomefile
open(15,file=nomefile,form='formatted')


do i=0,N-1

write (15,800) x(i),U(i)

end do
close(15)
write(*,*)'Fine scrittura su file ',nomefile
 

 
!----------------------------------------------
end do
!----------------------------------------------
close(16)


errore=maxval(Err)
err_rel=errore/2

write(*,*) 'N=.............................................',N
write(*,*) 'Errore=.....................................',errore
write(*,*) 'Errore relativo=............................',err_rel
write(3,500) errore,err_rel


500 format(2e16.7)
800 format (2e16.7)
close(3)






open(16,file='massimo.txt')

deallocate(x,U,Unext,F,Err,Xnew,vel)

allocate(xx_max(1:Ncicli),uu_max(1:Ncicli),vv_max(1:Ncicli),aa_max(1:Ncicli))

do i=1,Ncicli
read(16,*) tempo, xx_max(i), uu_max(i)
end do
close(16)
vv_max(1)=(xx_max(2)-xx_max(1))/(Nsalva*deltat)
aa_max(1)=(xx_max(3)-2*xx_max(2)+xx_max(1))/(Nsalva*deltat**2)
do i=2,Ncicli
vv_max(i)=(xx_max(i+1)-xx_max(i-1))/(2*Nsalva*deltat)
aa_max(i)=(xx_max(i+1)-2*xx_max(i)+xx_max(i-1))/(2*Nsalva*deltat**2)
end do
aa_max(Ncicli)=(xx_max(Ncicli)-2*xx_max(Ncicli-1)+xx_max(Ncicli-2))/(Nsalva*deltat**2)
vv_max(Ncicli)=(xx_max(Ncicli)-xx_max(Ncicli-1))/(Nsalva*deltat)

open(20,file='pos_val_max_vel_acc.txt')
do i=1,Ncicli

write(20,*) i*deltat, xx_max(i), uu_max(i), vv_max(i),aa_max(i)
end do
close(20)

deallocate(uu_max,vv_max,xx_max,aa_max)



 
end program kdvadim

!*********************************************************************
SUBROUTINE calcolaeffe(Us, Fs,N,deltaxs)
 integer :: i
 real, dimension(0:N-1) :: Us, Fs
 real :: deltaxs

 do i=2,N-2
   
  Fs(i)=(-1.5*(Us(i+1)**2-Us(i-1)**2)/deltaxs)-(Us(i+2)-2*Us(i+1)+2*Us(i-1)-Us(i-2))/(2*(deltaxs**3))
  
 end do
 
  Fs(N-1)=(-1.5*(Us(0)**2-Us(N-2)**2)/(deltaxs))-(Us(1)-2*Us(0)+2*Us(N-2)-Us(N-3))/(2*(deltaxs**3))
 
  Fs(N-2)=(-1.5*(Us(N-1)**2-Us(N-3)**2)/(deltaxs))-(Us(0)-2*Us(N-1)+2*Us(N-3)-Us(N-4))/(2*(deltaxs**3))

  Fs(0)=(-1.5*(Us(1)**2-Us(N-1)**2)/(deltaxs))-(Us(2)-2*Us(1)+2*Us(N-1)-Us(N-2))/(2*(deltaxs**3))  
  
  Fs(1)=(-1.5*(Us(2)**2-Us(0)**2)/(deltaxs))-(Us(3)-2*Us(2)+2*Us(0)-Us(N-1))/(2*(deltaxs**3))
  
END SUBROUTINE calcolaeffe
