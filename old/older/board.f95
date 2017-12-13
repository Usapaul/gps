program orbit
use newtonmethod
use ffunction
use ffunction1

implicit none

real(8), parameter :: mu=3.986004418d+14, c=299792458.0d0, OmegaE=7.2921151467d-5
integer, parameter :: typesobs=5 ! C1 P1 L1 D1 S1
integer, dimension(5) :: date, dateepoch
real(8) :: secdate, deltat, secepoch
integer :: i, j, numsat, numnavdata, id
integer, dimension(12) :: satnum
integer, dimension(40) :: idsave=0
type satellite
    sequence
    character(3) :: id
    real(8) :: C1, P1, L1, D1, S1
    real(8) :: poprav, hod, Vhoda
    real(8) :: IODE, Crs, Dn, M0
    real(8) :: Cuc, e, Cus, Asqrt
    real(8) :: epochsec, Cic, Omega0, Cis
    real(8) :: i0, Crc, omega, OmegaDOT
    real(8) :: IDOT, dop1, WN, dop2
    real(8) :: prec, state, Tgd, IODC
    real(8) :: tsend, dop3
    real(8) :: n0, P, tobs, tem, toe, t, n, Man, Ean, Van, argphi, du, dr, di
    real(8) :: u, r, iarg, Xorb, Yorb, Omegaist, Xzem, Yzem, Zzem
end type satellite
type(satellite), dimension(40) :: sat
real(8), dimension(4) :: XYZcdt0, XYZcdt
real(8) :: longitude, latitude, Teta, aell=6378137.0d0, bell=6356752.0d0, etrans, etrans1, rad, latdeg, longdeg
real(8), dimension(2,15) :: cords
real(8), dimension(2) :: cordsav
real(8), dimension(3) :: X1, X0
real(8) :: zhatie=298.257223563d0
rad=360.0d0/2.0d0/(4.0d0*atan(1.0d0))
bell=aell*(1.0d0-1.0d0/zhatie)
!---------------------------------------
sat%WN=0.0d0
forall (i=1:size(sat%id)) sat(i)%id='G'//achar(i/10+48)//achar(mod(i,10)+48)
open(100,file='obs.dat')
read(100,'(1x,i2.2,4(1x,i2),f11.7,2x,x,i3,12(x,i2),f12.9)') date, secdate, numsat, satnum, deltat
do i=1,numsat
    read(100,'(5(f14.3,x,x))') sat(satnum(i))%C1, sat(satnum(i))%P1, sat(satnum(i))%L1, sat(satnum(i))%D1, sat(satnum(i))%S1
enddo
close(100)
open(200,file='stationxxx.dat')
read(200,'(i2)') numnavdata
do i=1,numnavdata
    read(200,'(i2,1x,i2.2,4(1x,i2),f5.1,3d19.12)') id, dateepoch, secepoch, sat(id)%poprav, sat(id)%hod, sat(id)%Vhoda
    read(200,'(3x,4d19.12)') sat(id)%IODE, sat(id)%Crs, sat(id)%Dn, sat(id)%M0
    read(200,'(3x,4d19.12)') sat(id)%Cuc, sat(id)%e, sat(id)%Cus, sat(id)%Asqrt
    read(200,'(3x,4d19.12)') sat(id)%epochsec, sat(id)%Cic, sat(id)%Omega0, sat(id)%Cis
    read(200,'(3x,4d19.12)') sat(id)%i0, sat(id)%Crc, sat(id)%omega, sat(id)%OmegaDOT
    read(200,'(3x,4d19.12)') sat(id)%IDOT, sat(id)%dop1, sat(id)%WN, sat(id)%dop2
    read(200,'(3x,4d19.12)') sat(id)%prec, sat(id)%state, sat(id)%Tgd, sat(id)%IODC
    read(200,'(3x,4d19.12)') sat(id)%tsend, sat(id)%dop3
    idsave(i)=id
enddo
close(200)
!---------------------------------------
forall (i=1:size(sat%id)) sat(i)%toe=sat(i)%epochsec
forall (i=1:size(sat%id)) sat(i)%n0=sqrt(mu/sat(i)%Asqrt**6)
do i=1,size(sat%id)
if (sat(i)%WN>1.0d0) then
    sat(i)%tobs=mod((date(3)+30*(date(2)-6)-25)*24*3600+date(4)*3600+date(5)*60,604800)+secdate ! от 25 июня счет до июля
!    sat(i)%P=sat(i)%L1/(154.0d0*10.23d+6)*c
    sat(i)%P=sat(i)%P1
    call orbital(sat(i))
endif
enddo
forall (i=1:size(sat%id)) sat(i)%P=sat(i)%P+c*sat(i)%poprav!+sat(i)%hod*sat(i)%t)
!---------------------------------------
!go to 1000
!go to 1111
idsave(1:6)=(/2,4,12,25,29,31/)
do j=1,6
i=idsave(j)
write(*,*) 'name=', sat(i)%id
write(*,*) 'n0=', sat(i)%n0
write(*,*) 'tem=', sat(i)%tem
write(*,*) 't=', sat(i)%t
write(*,*) 'n=', sat(i)%n
write(*,*) 'M=', sat(i)%Man
write(*,*) 'E=', sat(i)%Ean
write(*,*) 'V=', sat(i)%Van
write(*,*) 'phi=', sat(i)%argphi
write(*,*) 'du=', sat(i)%du
write(*,*) 'dr=', sat(i)%dr
write(*,*) 'di=', sat(i)%di
write(*,*) 'u=', sat(i)%u
write(*,*) 'r=', sat(i)%r
write(*,*) 'i=', sat(i)%iarg
write(*,*) 'Xorb=', sat(i)%Xorb
write(*,*) 'Yorb=', sat(i)%Yorb
write(*,*) 'X=', sat(i)%Xzem
write(*,*) 'Y=', sat(i)%Yzem
write(*,*) 'Z=', sat(i)%Zzem
enddo
go to 1111

1000 i=12
write(*,*) 'id=', sat(i)%id
write(*,*) 'C1=', sat(i)%C1
write(*,*) 'P=', sat(i)%P
write(*,*) 'popravka=', sat(i)%poprav
write(*,*) 'tobs=', sat(i)%tobs
write(*,*) 'tem=', sat(i)%tem
write(*,*) 'toe=', sat(i)%toe
write(*,*) 't=', sat(i)%t
write(*,*) 'Asqrt=', sat(i)%Asqrt
write(*,*) 'n0=', sat(i)%n0
write(*,*) 'Dn=', sat(i)%Dn
write(*,*) 'n=', sat(i)%n
write(*,*) 'M0=', sat(i)%M0
write(*,*) 'Man=', sat(i)%Man
write(*,*) 'e=', sat(i)%e
write(*,*) 'Ean=', sat(i)%Ean
write(*,*) 'Van=', sat(i)%Van
write(*,*) 'omega=', sat(i)%omega
write(*,*) 'argphi=', sat(i)%argphi
write(*,*) 'Cus=', sat(i)%Cus, 'Cuc=', sat(i)%Cuc
write(*,*) 'Crs=', sat(i)%Crs, 'Crc=', sat(i)%Crc
write(*,*) 'Cis=', sat(i)%Cis, 'Cic=', sat(i)%Cic
write(*,*) 'du=', sat(i)%du
write(*,*) 'dr=', sat(i)%dr
write(*,*) 'di=', sat(i)%di
write(*,*) 'u=', sat(i)%u
write(*,*) 'r=', sat(i)%r
write(*,*) 'i0=', sat(i)%i0, 'IDOT=', sat(i)%IDOT
write(*,*) 'iarg=', sat(i)%iarg
write(*,*) 'Xorb=', sat(i)%Xorb
write(*,*) 'Yorb=', sat(i)%Yorb
write(*,*) 'Omega0=', sat(i)%Omega0
write(*,*) 'OmegaDOT=', sat(i)%OmegaDOT
write(*,*) 'Omegaist=', sat(i)%Omegaist
write(*,*) 'Xzem=', sat(i)%Xzem
write(*,*) 'Yzem=', sat(i)%Yzem
write(*,*) 'Zzem=', sat(i)%Zzem
!---------------------------------------
1111 open(800,file='saves.dat')
do i=1,15
read(800,*) idsave(1:4)
open(888,file='XYZPsat.dat')
do j=1,4
    write(888,*) sat(idsave(j))%Xzem
    write(888,*) sat(idsave(j))%Yzem
    write(888,*) sat(idsave(j))%Zzem
    write(888,*) sat(idsave(j))%P
enddo
close(888)
XYZcdt0=(/2765332.1835d0,1615597.8327d0,5497278.3855d0,c*deltat*1.0d-2/) ! начальное приближение
call newton(XYZcdt0,10000,XYZcdt,Ffun1)
etrans=sqrt((aell**2-bell**2)/aell**2)
etrans1=sqrt((aell**2-bell**2)/bell**2)
Teta=atan(aell/bell*XYZcdt(3)/sqrt(XYZcdt(1)**2+XYZcdt(2)**2))
longitude=atan(XYZcdt(2)/XYZcdt(1))
latitude=atan((XYZcdt(3)+etrans1**2*bell*(sin(Teta))**3)/(sqrt(XYZcdt(1)**2+XYZcdt(2)**2)-etrans**2*aell*(cos(Teta))**3))
latdeg=latitude*rad
longdeg=longitude*rad

cords(1,i)=latdeg
cords(2,i)=longdeg
!open(777,file='coords.dat',status='replace')
!write(777,*) latdeg
!write(777,*) longdeg
!close(777)
enddo
close(800)
!write(*,*) 'idsave=', idsave(1:6)
!write(*,'(6x,i2,x,i2,x,f13.10)') int(latdeg), int((latdeg-int(latdeg))*60), int(mod(latdeg*3600,60.0)*1000)*1.0/1000
!write(*,'(6x,i2,x,i2,x,f13.10)') int(longdeg), int((longdeg-int(longdeg))*60), int(mod(longdeg*3600,60.0)*1000)*1.0/1000
!write(*,'(6x,i2,x,i2,x,f13.10)') int(longdeg/15), int((longdeg/15-int(longdeg/15))*60), int(mod(longdeg/15*3600,60.0)*1000)*1.0/1000

cordsav(1)=sum(cords(1,1:14))/14
cordsav(2)=sum(cords(2,1:14))/14

write(*,*) (cordsav(1)-59.942214d0)*3600*30.92d0
write(*,*) (cordsav(2)-30.294900d0)*3600*30.92d0*0.5d0
write(*,*) '>>>'
do i=1,15
    write(*,*) (cords(1,i)-cordsav(1))*3600*30.92d0, (cords(2,i)-cordsav(2))*3600*30.92d0*0.5d0, cords(:,i)
enddo




!X0=XYZcdt0(1:3)
!do i=1,6
!    X1(1)=sat(idsave(i))%Xzem
!    X1(2)=sat(idsave(i))%Yzem
!    X1(3)=sat(idsave(i))%Zzem
!    write(*,*) rad*asin(sin(acos(dot_product(X1,X0)/dot_product(X1,X1)/dot_product(X0,X0)))*dot_product(X1,X1)*sqrt())
!enddo


contains

subroutine orbital(sat)
implicit none
type satellite
    sequence
    character(3) :: id
    real(8) :: C1, P1, L1, D1, S1
    real(8) :: poprav, hod, Vhoda
    real(8) :: IODE, Crs, Dn, M0
    real(8) :: Cuc, e, Cus, Asqrt
    real(8) :: epochsec, Cic, Omega0, Cis
    real(8) :: i0, Crc, omega, OmegaDOT
    real(8) :: IDOT, dop1, WN, dop2
    real(8) :: prec, state, Tgd, IODC
    real(8) :: tsend, dop3
    real(8) :: n0, P, tobs, tem, toe, t, n, Man, Ean, Van, argphi, du, dr, di
    real(8) :: u, r, iarg, Xorb, Yorb, Omegaist, Xzem, Yzem, Zzem
end type satellite
type(satellite) :: sat
character(3) :: id
real(8) :: C1, P1, L1, D1, S1
real(8) :: poprav, hod, Vhoda
real(8) :: IODE, Crs, Dn, M0
real(8) :: Cuc, e, Cus, Asqrt
real(8) :: epochsec, Cic, Omega0, Cis
real(8) :: i0, Crc, omega, OmegaDOT
real(8) :: IDOT, dop1, WN, dop2
real(8) :: prec, state, Tgd, IODC
real(8) :: tsend, dop3
real(8) :: n0, P, tobs, tem, toe, t, n, Man, Ean, Van, argphi, du, dr, di
real(8) :: u, r, iarg, Xorb, Yorb, Omegaist, Xzem, Yzem, Zzem

real(8), dimension(2) :: helpMan, helpEan
!---------------------------------------
C1=sat%C1; P1=sat%P1; L1=sat%L1; D1=sat%D1; S1=sat%S1
poprav=sat%poprav; hod=sat%hod; Vhoda=sat%Vhoda
IODE=sat%IODE; Crs=sat%Crs; Dn=sat%Dn; M0=sat%M0
Cuc=sat%Cuc; e=sat%e; Cus=sat%Cus; Asqrt=sat%Asqrt
epochsec=sat%epochsec; Cic=sat%Cic; Omega0=sat%Omega0; Cis=sat%Cis
i0=sat%i0; Crc=sat%Crc; omega=sat%omega; OmegaDOT=sat%OmegaDOT
IDOT=sat%IDOT; dop1=sat%dop1; WN=sat%WN; dop2=sat%dop2
prec=sat%prec; state=sat%state; Tgd=sat%Tgd; IODC=sat%IODC
tsend=sat%tsend; dop3=sat%dop3
n0=sat%n0; P=sat%P; tobs=sat%tobs; tem=sat%tem; toe=sat%toe; t=sat%t; n=sat%n
Man=sat%Man; Ean=sat%Ean; Van=sat%Van; argphi=sat%argphi; du=sat%du; dr=sat%dr; di=sat%di
u=sat%u; r=sat%r; iarg=sat%iarg; Xorb=sat%Xorb; Yorb=sat%Yorb
Omegaist=sat%Omegaist; Xzem=sat%Xzem; Yzem=sat%Yzem; Zzem=sat%Zzem
!---------------------------------------
tem=tobs-P/c
t=tem-toe
if (t>=302400.0d0) then
    t=t-604800
elseif (t<=-302400.0d0) then
    t=t+604800
endif
n=n0+Dn
Man=M0+n*t
helpMan=Man
open(555,file='eMan.dat')
write(555,*) e
write(555,*) Man
close(555)
call newton(helpMan,100000,helpEan,Ffun)
Ean=helpEan(1)
Van=atan(sqrt((1+e)/(1-e))*tan(Ean/2))*2
argphi=Van+omega
du=Cus*sin(2*argphi)+Cuc*cos(2*argphi)
dr=Crs*sin(2*argphi)+Crc*cos(2*argphi)
di=Cis*sin(2*argphi)+Cic*cos(2*argphi)
u=argphi+du
r=Asqrt**2*(1.0d0-e*cos(Ean))+dr
iarg=i0+di+IDOT*t
Xorb=r*cos(u)
Yorb=r*sin(u)
Omegaist=Omega0+(OmegaDOT-OmegaE)*t-OmegaE*toe
Xzem=Xorb*cos(Omegaist)-Yorb*cos(iarg)*sin(Omegaist)
Yzem=Xorb*sin(Omegaist)+Yorb*cos(iarg)*cos(Omegaist)
Zzem=Yorb*sin(iarg)
!---------------------------------------
sat%P=P; sat%tobs=tobs; sat%tem=tem; sat%toe=toe; sat%t=t; sat%n=n
sat%Man=Man; sat%Ean=Ean; sat%Van=Van; sat%argphi=argphi; sat%du=du; sat%dr=dr; sat%di=di
sat%u=u; sat%r=r; sat%iarg=iarg; sat%Xorb=Xorb; sat%Yorb=Yorb
sat%Omegaist=Omegaist; sat%Xzem=Xzem; sat%Yzem=Yzem; sat%Zzem=Zzem



end subroutine orbital



end program orbit