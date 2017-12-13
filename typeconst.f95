module typeconst

implicit none

integer, parameter  :: pr=8 ! Константа, задающая разновидность real

	! Значения параметров вращения Земли и скорости света
real(pr), parameter :: mu=3.986004418d+14, c=299792458.0d0, OmegaE=7.2921151467d-5
	! Количество представленных...
	! > наблюдений => Nobs
	! > типов наблюдений => Ntype
	! > спутников с данными => Nsat
	! И количество строк в заголовке rinex файла
integer :: Nobs, Ntype, Nsat, Nheaderlines
	! В заголовке указаны приблизительные координаты. Они записываются в этот массив 
real(pr), dimension(3) :: approx_XYZ


! Нужные данные:
! > Дата и время наблюдений -- одномерный массив из переменных типа
! type(datewithsec) с размером = число наблюдений Nobs
! > Названия типов наблюдений -- одномерный массив из переменных типа
! character(2) с размером = число типов Ntype
! > Имена спутников -- одномерный массив из переменных типа 
! character(3) с размером = число спутников Nsat



type datewithsec
	integer, dimension(5) :: datehm
	real(pr) :: sec
end type datewithsec

type obs_data
	real(pr) :: value
	integer  :: LLI, power
end type obs_data

type obs_firstline
	type(datewithsec) :: datetime
	integer :: flag, numsat
	real(pr) :: delta_t
end type obs_firstline

!type one_sat
!	character(1) :: idx
!	character(2) :: idn
!	
!end type one_sat

end module typeconst