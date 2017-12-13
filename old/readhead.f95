program readhead

implicit none

integer :: i, j, k, nheaderlines
character(20), dimension(0:500) :: headermarks ! Предполагается, что в заголовке не больше 500 строк
real(8) :: FormatOfVersion
character(1) :: TypeOfFile, SatelliteSystem
character(20) :: NameOfProgram, NameOfAgency, DateOfCreation
character(60) :: identcomment
character(60) :: MarkerName
character(20) :: MarkerNumber
character(20) :: ObserverName
character(40) :: AgencyName
character(20) :: NumberOfRec, TypeOfRec, VersionOfRec
character(20) :: NumberOfAntenna, TypeOfAntenna
real(8), dimension(3) :: ApproxCoordinates, DeltaHEN
integer, dimension(3) :: L1L2and0
!--------------------------------------
integer :: L1second, L2second, NumberOfGoodSatellite
integer :: helpL1L2=0
type PRN
    sequence
    character(1) :: idchar
    integer :: idnum
    character(2), dimension(27) :: TypeObs
    integer, dimension(27) :: NumObs
end type PRN
type(PRN), dimension(7) :: ListPRN
!--------------------------------------
integer :: helpTypesObs=0, TypesObsNum
type ObsType
    sequence
    character(1) :: iddata
    character(1) :: idfreq
end type ObsType
type(ObsType), dimension(27) :: ObsTypes
!--------------------------------------
real(8) :: IntervalObs
integer, dimension(5) :: TimeFirst, TimeLast
real(8) :: TimeFirstSec, TimeLastSec
character(3) :: SystemOfTime
integer :: RCVClock, LeapSeconds, NumOfSat
type(PRN), allocatable :: SatTotal(:)
integer :: helpPRN=0

!---------------------------------------

open(200,file='rinexdata.dat') ! Сначала производится чтение в файле только меток в заголовке
i=0
do while ( headermarks(i) .ne. 'END OF HEADER       ' )
    i=i+1
    read(200,'(60x,a20)') headermarks(i)
enddo
nheaderlines=i
close(200)

open(100,file='rinexdata.dat') ! По известным меткам производится считывание данных из конкретной строки
do k=1,nheaderlines
    select case(headermarks(k))
    case('RINEX VERSION / TYPE')
    read(100,'(f9.2,11x,a1,19x,a1,19x)') FormatOfVersion, TypeOfFile, SatelliteSystem
    case('PGM / RUN BY / DATE ')
    read(100,'(a20,a20,a20)') NameOfProgram, NameOfAgency, DateOfCreation
    case('COMMENT             ')
    read(100,'(a60)') identcomment
    case('MARKER NAME         ')
    read(100,'(a20)') MarkerName
    case('MARKER NUMBER       ')
    read(100,'(a20)') MarkerNumber
    case('OBSERVER / AGENCY   ')
    read(100,'(a20,a40)') ObserverName, AgencyName
    case('REC # / TYPE / VERS ')
    read(100,'(3(a20))') NumberOfRec, TypeOfRec, VersionOfRec
    case('ANT # / TYPE        ')
    read(100,'(2(a20))') NumberOfAntenna, TypeOfAntenna
    case('APPROX POSITION XYZ ')
    read(100,'(3(f14.4))') ApproxCoordinates(1:3)
    case('ANTENNA: DELTA H/E/N')
    read(100,'(3(f14.4))') DeltaHEN(1:3)
    case('WAVELENGTH FACT L1/2')
    if (helpL1L2==0) then                    ! helpL1L2 используется только потому,
        read(100,'(3(i6))') L1L2and0(1:3)    ! что в rinex файле присутствуют две
        helpL1L2=1                           ! одинаковые строки WAVELENGTH FACT,
    else                                     ! и данные в них считываются по очереди
        read(100,'(3(i6),7(3x,a1,i2))') L1second, L2second, NumberOfGoodSatellite, (ListPRN(i)%idchar,ListPRN(i)%idnum,i=1,7)
    endif
    case('# / TYPES OF OBSERV ')
    ObsTypes%iddata=' '
    ObsTypes%idfreq=' '
    if (helpTypesObs==0) then
    ! Типов может быть больше 9, условно я выбрал, что уж точно меньше 27.
    ! Если больше 9, то запись производится в следующие строки.
    ! Учет возможного наличия нескольких строк с # / TYPES OF OBSERV происходит посредством переменной helpTypesObs
    ! То есть, если select case попадает сюда снова, то считываются данные со второй (третьей) строки.
        read(100,'(i6,9(4x,a1,a1))') TypesObsNum, (ObsTypes(i)%iddata,ObsTypes(i)%idfreq,i=1,9)
        helpTypesObs=helpTypesObs+1
    elseif (helpTypesObs==1) then
        read(100,'(6x,9(4x,a1,a1))') (ObsTypes(i)%iddata,ObsTypes(i)%idfreq,i=10,18)
        helpTypesObs=helpTypesObs+1
    else
        read(100,'(6x,9(4x,a1,a1))') (ObsTypes(i)%iddata,ObsTypes(i)%idfreq,i=19,27)
        helpTypesObs=helpTypesObs+1
    endif
    helpTypesObs=helpTypesObs-1
    case('INTERVAL            ')
    read(100,'(f10.3)') IntervalObs
    case('TIME OF FIRST OBS   ')
    read(100,'(5i6,f13.7,5x,a3)') TimeFirst(1:5),TimeFirstSec, SystemOfTime
    case('TIME OF LAST OBS    ')
    read(100,'(5i6,f13.7,5x,a3)') TimeLast(1:5), TimeLastSec, SystemOfTime
    case('RCV CLOCK OFFS APPL ')
    read(100,'(i6)') RCVClock
    case('LEAP SECONDS        ')
    read(100,'(i6)') LeapSeconds
    case('# OF SATELLITES     ')
    read(100,'(i6)') NumOfSat
    allocate (SatTotal(NumOfSat))
    case('PRN / # OF OBS      ')
    if (helpPRN==0) then
    forall (i=1:27) SatTotal%TypeObs(i)=ObsTypes(i)%iddata//ObsTypes(i)%idfreq
    do j=1,NumOfSat
        read(100,'(3x,a1,i2,9i6)') SatTotal(j)%idchar, SatTotal(j)%idnum, SatTotal(j)%NumObs(1:9)
        if (helpTypesObs>0) read(100,'(6x,9i6)') SatTotal(j)%NumObs(10:18)
        if (helpTypesObs>1) read(100,'(6x,9i6)') SatTotal(j)%NumObs(19:27)
    enddo
    helpPRN=helpPRN+1 ! Просто +1, чтобы условие if не выполнялось, так как за один раз считываются все строки PRN / # OF OBS
    endif
    case('END OF HEADER       ')
    continue
    case default
    write(*,*) 'Unknown string:', headermarks(k)
    end select
enddo
close(100)
!--------------------------------------
open(333,file='varsavemod.f95')
    write(333,*) 'module varsavemod'
    write(333,*) 'integer, parameter :: NumOfSat=',     NumOfSat
    write(333,*) 'integer, parameter :: TypesObsNum=',  TypesObsNum
    write(333,*) 'integer, parameter :: nheaderlines=', nheaderlines
    write(333,*) 'real(8), parameter :: X00=', ApproxCoordinates(1)
    write(333,*) 'real(8), parameter :: Y00=', ApproxCoordinates(2)
    write(333,*) 'real(8), parameter :: Z00=', ApproxCoordinates(3)   
    write(333,*) 'end module varsavemod'
close(333)

open(900,file='namestypes.dat')
    do i=1,TypesObsNum
        write(900,'(a2)') SatTotal(1)%TypeObs(i)
    enddo
close(900)

open(400,file='namessat.dat')
    do j=1,NumOfSat
        write(400,'(a1,i2)') SatTotal(j)%idchar, SatTotal(j)%idnum
    enddo
close(400)

end program readhead







