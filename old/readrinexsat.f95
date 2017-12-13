program readfile
use varsavemod

use newtonmethod
use ffunction
use ffunction1

implicit none

type observation
    sequence
    integer, dimension(5) :: date ! В эти массив и переменную
    real(8)               :: sec  ! записывается дата и время наблюдения
    integer :: flag
    integer :: numsat              
    character(1), dimension(((NumOfSat-1)/12+1)*12) :: idx ! Всего спутников с данными наблюдений указано в заголовке и считывается в readhead.
    integer,      dimension(((NumOfSat-1)/12+1)*12) :: idn ! Это число равно NumOfSat, но для чтения rinex лучше число, кратное 12, большее NumOfSat 
    real(8) :: dt
    real(8), dimension(30,((NumOfSat-1)/12+1)*12) :: obsdata ! Для каждого из NumOfSat спутников в строке матрицы -- данные наблюдений
    integer, dimension(30,((NumOfSat-1)/12+1)*12) :: LLI     ! 30 > 27 и кратно 5, это тоже нужно для удобного считывания
    integer, dimension(30,((NumOfSat-1)/12+1)*12) :: power
end type observation
type(observation), allocatable :: list(:)
!--------------------------------------
character(2), dimension(TypesObsNum) :: typesnames
!--------------------------------------
integer :: xtype=1 ! будет выбран конкретный тип наблюдений получения псевдодальности
integer :: xobs=1 ! будет выбран конкретный номер наблюдения
!--------------------------------------
integer :: Ntypes=TypesObsNum, Nobs
integer :: i, j, k, n, m, nlinesfor1sat
character(3)  :: jx
character(1)  :: dump
character(20) :: dumpmed
character(80) :: dumpbig
! Предполагается, что наблюдений не больше 10000. В массив NumbersOfSat 
! записывается количество спутников с данными для каждого момента наблюдений
integer, dimension(10000) :: NumbersOfSat=0
!--------------------------------------
type satellite
    sequence
    character(1), dimension(NumOfSat) :: idx
    integer, dimension(NumOfSat)      :: idn
    integer, allocatable :: date(:,:)
    real(8), allocatable :: sec(:)
    real(8), allocatable :: data(:,:,:)
end type satellite
type(satellite) :: sat

!--------------------------------------

nlinesfor1sat=(Ntypes-1)/5+1 ! nlinesfor1sat обозначает количество строк, отведенных под запись наблюдений для одного спутника

open(101,file='rinexdata.dat')
    dump=' '
    do i=1,nheaderlines
        read(101,'(a1)') dump
    enddo
    Nobs=0 ! 0 -- это просто стартовое число, оно будет увеличиваться по ходу чтения
    do k=1,size(NumbersOfSat)
        Nobs=Nobs+1
        read(101,'(29x,i3)',end=88) NumbersOfSat(Nobs)
        do j=2,(NumbersOfSat(Nobs)-1)/12+1
            read(101,'(a1)') dump
        enddo
        do i=1,NumbersOfSat(Nobs)
            do j=1,nlinesfor1sat
                read(101,'(a1)') dump
            enddo
        enddo
    enddo
    88 Nobs=Nobs-1
close(101)
!--------------------------------------
allocate (list(1:Nobs))

open(333,file='rinexdata.dat')
    do i=1,nheaderlines
        read(333,'(a1)') dump
    enddo
    do k=1,Nobs
        read(333,1111) list(k)%date,list(k)%sec,list(k)%flag,list(k)%numsat,(list(k)%idx(j),list(k)%idn(j),j=1,12),list(k)%dt
        1111 format (1x,i2.2,4(1x,i2),f11.7,2x,i1,i3,12(a1,i2),f12.9)
        do i=2,(list(k)%numsat-1)/12+1
            read(333,'(32x,12(a1,i2))') (list(k)%idx(j),list(k)%idn(j),j=(i-1)*12+1,i*12) 
        enddo
        n=NumbersOfSat(k)    
        do j=1,n
            do i=1,nlinesfor1sat              
                read(333,'(5(f14.3,i1,i1))') (list(k)%obsdata(m,j),list(k)%LLI(m,j),list(k)%power(m,j),m=(i-1)*5+1,i*5)
            enddo
        enddo
    enddo
close(333)
!--------------------------------------
!allocate (sat%date(5,Nobs,NumOfSat),sat%sec(Nobs,NumOfSat),sat%data(30,Nobs,NumOfSat))
allocate (sat%data(30,Nobs,NumOfSat))
allocate (sat%date(5,Nobs),sat%sec(Nobs))
forall (k=1:Nobs)
    sat%date(1:5,k)=list(k)%date(1:5)
    sat%sec(k)=list(k)%sec
end forall

open(500,file='namessat.dat')
    do j=1,NumOfSat
        read(500,'(a1,i2)') sat%idx(j), sat%idn(j)
    enddo
close(500)

do j=1,NumOfSat
    do k=1,Nobs
        do i=1,list(k)%numsat
            if ( (sat%idx(j) .eq. list(k)%idx(i)) .and. (sat%idn(j) .eq. list(k)%idn(i)) ) then
                sat%data(1:30,k,j)=list(k)%obsdata(1:30,i)
            endif
        enddo
    enddo
enddo


!xtype=1 ! здесь можно в цикле запускать процедуру для разных типов
!xobs=1  ! и разных порядковых номерах наблюдений
call findsatandme(xtype,xobs)

!--------------------------------------
! В следующих блоках производится чтение навигационного файла и расчет координат
!--------------------------------------
contains

subroutine findsatandme(xt,xo)

implicit none

real(8), parameter :: mu=3.986004418d+14, c=299792458.0d0, OmegaE=7.2921151467d-5
type getorbit
    sequence
    character(1) :: idx
    integer      :: idn
    integer, dimension(5) :: epochdate
    real(8)               :: epochsec
    real(8) :: poprav, Vhod, axhoda
    real(8) :: IODE, Crs, Dn, M0
    real(8) :: Cuc, e, Cus, Asqrt
    real(8) :: toe0, Cic, Omega0, Cis
    real(8) :: i0, Crc, omega, OmegaDOT
    real(8) :: IDOT, dop1, WN, dop2
    real(8) :: prec, state, Tgd, IODC
    real(8) :: tsend, dop3
    real(8) :: n0, P, tobs, tem, toe, t, n
    real(8) :: Man, Ean, Van, argphi, du, dr, di
    real(8) :: u, r, iarg, Xorb, Yorb
    real(8) :: Omegaist, Xzem, Yzem, Zzem
end type getorbit
type(getorbit), dimension(100) :: orbit ! Условно предполагаем, что номера спутников (GPS) не превосходят 100    
!--------------------------------------
integer :: i, j, k, n, m, nexist=0 ! nexist -- число спутников, для которых есть навигационные данные
integer, dimension(NumOfSat) :: idsave
integer, intent(in) :: xt ! тип наблюдений для расчета псевдодальности
integer, intent(in) :: xo ! порядковый номер наблюдения
integer :: week ! Текущая неделя GPS момента наблюдений, понадобится в цикле
real(8) :: nsec
integer :: Cnk ! Коэффициент для расчета количества перестановок
integer, allocatable :: combs(:,:)
real(8) :: latitude, longitude
real(8), allocatable :: cords(:,:)
!--------------------------------------

!forall (j=1:NumOfSat) orbit(j)%idx=sat%idx(j)
!forall (j=1:NumOfSat) orbit(j)%idn=sat%idn(j)

forall (j=1:size(orbit%idx)) orbit(j)%idx='G'
forall (j=1:size(orbit%idx)) orbit(j)%idn=j
orbit%WN=0.0d0

open(303,file='gpsnavdata.dat')
    dumpmed=' '
    do while (dumpmed .ne. 'END OF HEADER       ')
        read(303,'(60x,a20)') dumpmed
    enddo
    do i=1,NumOfSat
        read(303,2222,end=99) j, orbit(j)%epochdate, orbit(j)%epochsec, orbit(j)%poprav, orbit(j)%Vhod, orbit(j)%axhoda
        2222 format (i2,x,i2.2,4(1x,i2),f5.1,3d19.12) 
        read(303,3333) orbit(j)%IODE,  orbit(j)%Crs,   orbit(j)%Dn,     orbit(j)%M0
        read(303,3333) orbit(j)%Cuc,   orbit(j)%e,     orbit(j)%Cus,    orbit(j)%Asqrt
        read(303,3333) orbit(j)%toe0,  orbit(j)%Cic,   orbit(j)%Omega0, orbit(j)%Cis
        read(303,3333) orbit(j)%i0,    orbit(j)%Crc,   orbit(j)%omega,  orbit(j)%OmegaDOT
        read(303,3333) orbit(j)%IDOT,  orbit(j)%dop1,  orbit(j)%WN,     orbit(j)%dop2
        read(303,3333) orbit(j)%prec,  orbit(j)%state, orbit(j)%Tgd,    orbit(j)%IODC
        read(303,3333) orbit(j)%tsend, orbit(j)%dop3
        3333 format (3x,4d19.12)
        nexist=nexist+1
        idsave(i)=j
    enddo
    99 continue
close(303)
!--------------------------------------
forall (i=1:size(orbit%idx)) orbit(i)%toe=orbit(i)%toe0
forall (i=1:size(orbit%idx)) orbit(i)%n0=sqrt(mu/orbit(i)%Asqrt**6)

do i=1,size(orbit%idx)
    !call numWN(sat%date(1:3,xo),week)
    !if (abs(orbit(i)%WN-week)>1) then
    !    write(*,*) i, orbit(i)%WN, week
    !    write(*,*) 'Different weeks i rinex & NAVDATA. Continue? y/n'
    !    read(*,*) dump
    !    if (dump .eq. 'y') then
    !        continue
    !    else
    !        stop 'WN!=obs_week; Check your data and try again.'
    !    endif
    !endif
    !call secbtw(sat%date(1:5,xo),sat%sec(xo),orbit(i)%epochdate(1:5),orbit(i)%epochsec,nsec)
    !orbit(i)%tobs=nsec
    !call shortnames(orbit(i),xt,xo)
    if (orbit(i)%WN>1.0d0) then
        call secbtw(sat%date(1:5,xo),sat%sec(xo),orbit(i)%epochdate(1:5),orbit(i)%epochsec,nsec)
        orbit(i)%tobs=nsec+orbit(i)%toe0
        call shortnames(orbit(i),xt,xo)
    endif
enddo
forall (i=1:size(orbit%idx)) orbit(i)%P=orbit(i)%P+c*orbit(i)%poprav
!--------------------------------------
Cnk=product((/(i,i=1,nexist)/))/(1*2*3*4)/product((/(i,i=1,nexist-4)/)) ! Здесь записан факториал
allocate (combs(4,Cnk))
select case (nexist)
    case(1:3)
        stop 'The number of satellites with NAVDATA less than 4'
    case(0)
        stop 'No satellites with NAVDATA'
    case(4:)
        call getcombinations(nexist,Cnk,combs)
end select


allocate (cords(2,Cnk))

do k=1,Cnk
    open(900,file='XYZPsat.dat')
    do j=1,4
        write(900,*) orbit(idsave(combs(j,k)))%Xzem
        write(900,*) orbit(idsave(combs(j,k)))%Yzem
        write(900,*) orbit(idsave(combs(j,k)))%Zzem
        write(900,*) orbit(idsave(combs(j,k)))%P
    enddo
    close(900)
    call getXYZ(latitude,longitude)
    cords(1,k)=latitude
    cords(2,k)=longitude
enddo

write(jx,'(i1,i2)') xt,xo

open(707,file='cords'//jx//'.dat')
    do k=1,Cnk
        write(707,*) cords(1:2,k)
    enddo
close(707)


end subroutine findsatandme


subroutine numWN(date,nw)
! В процедуре производится вычисление номера недели nw на текущую дату date(3).date(2).date(1)

implicit none

integer, dimension(3), parameter   :: week0=(/80, 1, 6/)
integer, parameter                 :: to00=7300 ! 7300 дней с 06.01.1980 по 01.01.2000
integer, dimension(3), intent(in)  :: date
integer,               intent(out) :: nw
!--------------------------------------
integer :: i, days
integer, dimension(12) :: dinmonths=(/31,28,31,30,31,30,31,31,30,31,30,31/)

!--------------------------------------

days=to00 ! сразу переходим к 21 веку
days=days-1 ! чтобы все-таки вернуться именно на 31.12.1999, так удобнее считать

do i=1,date(1)  
    days=days+365
    if (mod(date(1),4)==0) days=days+1
enddo    


do i=1,date(2)-1
    days=days+dinmonths(i)
    if (i==3 .and. mod(date(1),4)==0) days=days+1 ! учтем, что в високосном году в феврале 29 дней
enddo    

days=days+date(3) ! Наконец, пришли к нужной дате
days=days-1 ! Потому что считали количество прошедших с 06.01.1980 дней, а нужен номер недели
            ! То есть прошло 7 дней, вроде бы неделя, но еще + время, и уже, значит, следующая неделя
nw=(days-1)/7+1

end subroutine numWN


subroutine secbtw(date2,sec2,date1,sec1,secs)

implicit none

integer, dimension(5), intent(in)  :: date1, date2
real(8),               intent(in)  :: sec1, sec2               
real(8),               intent(out) :: secs
integer, parameter                 :: to00=7300-1 ! 7300 дней с 06.01.1980 по 01.01.2000
!--------------------------------------
integer :: i, days1, days2
integer, dimension(12) :: dinmonths=(/31,28,31,30,31,30,31,31,30,31,30,31/)

!--------------------------------------

days1=to00
days2=to00

do i=1,date1(1)  
    days1=days1+365
    if (mod(date1(1),4)==0) days1=days1+1
enddo 
do i=1,date2(1)  
    days2=days2+365
    if (mod(date2(1),4)==0) days2=days2+1
enddo 

do i=1,date1(2)-1
    days1=days1+dinmonths(i)
    if (i==3 .and. mod(date1(1),4)==0) days1=days1+1 ! учтем, что в високосном году в феврале 29 дней
enddo    
do i=1,date2(2)-1
    days2=days2+dinmonths(i)
    if (i==3 .and. mod(date2(1),4)==0) days2=days2+1 ! учтем, что в високосном году в феврале 29 дней
enddo    

days1=days1+date1(3) ! Наконец, пришли к нужной дате
!days1=days1-1 ! Потому что считали количество прошедших с 06.01.1980 дней, а нужен номер недели
            ! То есть прошло 7 дней, вроде бы неделя, но еще + время, и уже, значит, следующая неделя
days2=days2+date2(3) ! Наконец, пришли к нужной дате
!days2=days2-1 

secs=(days2-days1)*86400+(date2(4)-date1(4))*3600+(date2(5)-date1(5))*60+sec2-sec1

end subroutine secbtw

subroutine shortnames(orbit,xt,xo)

implicit none

integer, intent(in) :: xt, xo
real(8), parameter :: mu=3.986004418d+14, c=299792458.0d0, OmegaE=7.2921151467d-5
!--------------------------------------
type getorbit
    sequence
    character(1) :: idx
    integer      :: idn
    integer, dimension(5) :: epochdate
    real(8)               :: epochsec
    real(8) :: poprav, Vhod, axhoda
    real(8) :: IODE, Crs, Dn, M0
    real(8) :: Cuc, e, Cus, Asqrt
    real(8) :: toe0, Cic, Omega0, Cis
    real(8) :: i0, Crc, omega, OmegaDOT
    real(8) :: IDOT, dop1, WN, dop2
    real(8) :: prec, state, Tgd, IODC
    real(8) :: tsend, dop3
    real(8) :: n0, P, tobs, tem, toe, t, n
    real(8) :: Man, Ean, Van, argphi, du, dr, di
    real(8) :: u, r, iarg, Xorb, Yorb
    real(8) :: Omegaist, Xzem, Yzem, Zzem
end type getorbit
type(getorbit) :: orbit    
!--------------------------------------
real(8) :: poprav, Vhod, axhoda
real(8) :: IODE,   Crs,   Dn,     M0
real(8) :: Cuc,    e,     Cus,    Asqrt
real(8) :: toe0,   Cic,   Omega0, Cis
real(8) :: i0,     Crc,   omega,  OmegaDOT
real(8) :: IDOT,   dop1,  WN,     dop2
real(8) :: prec,   state, Tgd,    IODC
real(8) :: tsend,  dop3

real(8) :: n0,       P,    tobs, tem,    toe, t,  n        
real(8) :: Man,      Ean,  Van,  argphi, du,  dr, di       
real(8) :: u,        r,    iarg, Xorb,   Yorb           
real(8) :: Omegaist, Xzem, Yzem, Zzem                   
!--------------------------------------
integer :: i, j
real(8), dimension(2) :: helpMan, helpEan

!--------------------------------------
poprav=orbit%poprav; Vhod=orbit%Vhod  ; axhoda=orbit%axhoda                                                                                 
IODE  =orbit%IODE  ; Crs  =orbit%Crs  ; Dn    =orbit%Dn    ; M0      =orbit%M0                                                                       
Cuc   =orbit%Cuc   ; e    =orbit%e    ; Cus   =orbit%Cus   ; Asqrt   =orbit%Asqrt                                                                      
toe0  =orbit%toe0  ; Cic  =orbit%Cic  ; Omega0=orbit%Omega0; Cis     =orbit%Cis                                                                        
i0    =orbit%i0    ; Crc  =orbit%Crc  ; omega =orbit%omega ; OmegaDOT=orbit%OmegaDOT                                                                  
IDOT  =orbit%IDOT  ; dop1 =orbit%dop1 ; WN    =orbit%WN    ; dop2    =orbit%dop2                                                                    
prec  =orbit%prec  ; state=orbit%state; Tgd   =orbit%Tgd   ; IODC    =orbit%IODC                                                                      
tsend =orbit%tsend ; dop3 =orbit%dop3  

n0      =orbit%n0      ; P   =orbit%P   ; tobs=orbit%tobs; tem   =orbit%tem   ; toe =orbit%toe ; t =orbit%t ; n =orbit%n   
Man     =orbit%Man     ; Ean =orbit%Ean ; Van =orbit%Van ; argphi=orbit%argphi; du  =orbit%du  ; dr=orbit%dr; di=orbit%di   
u       =orbit%u       ; r   =orbit%r   ; iarg=orbit%iarg; Xorb  =orbit%Xorb  ; Yorb=orbit%Yorb         
Omegaist=orbit%Omegaist; Xzem=orbit%Xzem; Yzem=orbit%Yzem; Zzem  =orbit%Zzem    
!--------------------------------------
do j=1,NumOfSat
    if ( (orbit%idx .eq. sat%idx(j)) .and. (orbit%idn .eq. sat%idn(j)) ) then
        P=sat%data(xt,xo,j)
    endif
enddo

tem=tobs-P/c
t=tem-toe
n=n0+Dn
Man=M0+n*t
helpMan=Man

open(550,file='eMan.dat')
    write(550,*) e
    write(550,*) Man
close(550)

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
!--------------------------------------
orbit%n0      =n0      ; orbit%P   =P   ; orbit%tobs=tobs; orbit%tem   =tem   ; orbit%toe =toe ; orbit%t =t ; orbit%n =n   
orbit%Man     =Man     ; orbit%Ean =Ean ; orbit%Van =Van ; orbit%argphi=argphi; orbit%du  =du  ; orbit%dr=dr; orbit%di=di   
orbit%u       =u       ; orbit%r   =r   ; orbit%iarg=iarg; orbit%Xorb  =Xorb  ; orbit%Yorb=Yorb         
orbit%Omegaist=Omegaist; orbit%Xzem=Xzem; orbit%Yzem=Yzem; orbit%Zzem  =Zzem             

end subroutine shortnames


subroutine getcombinations(n,Cnk,u)

implicit none

integer,                   intent(in)  :: n, Cnk
integer, dimension(4,Cnk), intent(out) :: u
!--------------------------------------
integer, allocatable :: s(:,:)
integer :: i, j, k, m, nbig
integer :: i0, i1, i2, i3
!--------------------------------------

nbig=product((/(n-i,i=0,3)/))
allocate (s(4,nbig))
!Cnk=product((/(i,i=1,n)/))/(1*2*3*4)/product((/(i,i=1,n-4)/))
!allocate (u(4,Cnk))
!--------------------------------------
k=0
do i0=1,n
    do i1=1,n 
        do i2=1,n
            do i3=1,n
                if ( (i3 .ne. i2) .and. (i3 .ne. i1) .and. (i3 .ne. i0) ) then
                    if ( (i2 .ne. i1) .and. (i2 .ne. i0) ) then
                        if ( (i1 .ne. i0) ) then
                            k=k+1
                            s(1:4,k)=(/i0,i1,i2,i3/)
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
enddo
!--------------------------------------
j=0
do k=1,nbig
    j=j+1
    do i=1,k-1
        do i0=1,4
            if (s(1,k)==s(i0,i)) then
                do i1=1,4
                    if (s(2,k)==s(i1,i)) then
                        do i2=1,4
                            if (s(3,k)==s(i2,i)) then
                                do i3=1,4
                                    if (s(4,k)==s(i3,i)) then
                                        s(:,k)=0
                                        go to 55
                                    endif
                                enddo
                            endif
                        enddo
                    endif
                enddo
            endif
        enddo
    enddo
    go to 77
    55 j=j-1
    go to 99
    77 u(:,j)=s(:,k)
    99 continue
enddo

!do i=1,Cnk
!    write(*,'(2x,4(i2,x))') u(:,i)
!enddo

end subroutine getcombinations


subroutine getXYZ(lat,long)

implicit none

real(8), parameter :: mu=3.986004418d+14, c=299792458.0d0, OmegaE=7.2921151467d-5
real(8), parameter ::  aell=6378137.0d0, bell=6356752.0d0, rad=360.0d0/2.0d0/(4.0d0*atan(1.0d0))
real(8), dimension(4) :: XYZcdt0, XYZcdt
real(8) :: etrans, etrans1


real(8) :: long, lat, Teta, latdeg, longdeg
real(8), dimension(2,15) :: cords
!real(8), dimension(2) :: cordsav
!real(8), dimension(3) :: X1, X0

!--------------------------------------

XYZcdt0(1)=X00
XYZcdt0(2)=Y00
XYZcdt0(3)=Z00
XYZcdt0(4)=list(xobs)%dt*c*1.0d-2
call newton(XYZcdt0,10000,XYZcdt,Ffun1)
!--------------------------------------
etrans=sqrt((aell**2-bell**2)/aell**2)
etrans1=sqrt((aell**2-bell**2)/bell**2)
Teta=atan(aell/bell*XYZcdt(3)/sqrt(XYZcdt(1)**2+XYZcdt(2)**2))
long=atan(XYZcdt(2)/XYZcdt(1))
lat=atan((XYZcdt(3)+etrans1**2*bell*(sin(Teta))**3)/(sqrt(XYZcdt(1)**2+XYZcdt(2)**2)-etrans**2*aell*(cos(Teta))**3))
lat=lat*rad
long=long*rad

end subroutine getXYZ

!--------------------------------------
!open(205,file='namestypes.dat')
!    do i=1,Ntypes
!        read(205,'(a2)') typesnames(i)
!    enddo
!close(205)
!
!do j=1,NumOfSat
!    write(jx,'(a1,i2)') sat%idx(j), sat%idn(j)
!    open(800,file=jx//'.dat')
!        write(800,'("Date & time of observation",x,*(8x,a2,8x))') typesnames
!        do k=1,Nobs
!            do i=1,list(k)%numsat
!                if ( (sat%idx(j) .eq. list(k)%idx(i)) .and. (sat%idn(j) .eq. list(k)%idn(i)) ) then
!                    write(800,'(x,i2.2,4(x,i2),f11.7,*(2x,f14.3,2x))') list(k)%date, list(k)%sec, list(k)%obsdata(1:Ntypes,i)
!                endif
!            enddo
!        enddo
!    close(800)
!enddo


end program readfile