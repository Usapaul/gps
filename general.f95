program getall
use typeconst
use readingheadermod
use readingobsmod

implicit none
integer, parameter :: nout=5
integer     :: i, j, k, s
integer, dimension(:,:), allocatable :: megaprov, justsats
real(pr) :: eps=1.e-3
type(obs_data),    dimension(:,:,:), allocatable :: datavalues0

call readingheader()

if (allocated(sat_names)) then          ! 
	call readingobs()                   ! Если в заголовке не было указано количество спутников с данными, то    
else if (Nsat>0) then                   ! процедура readingobs будет вызвана с учетом этого. Если было указано количество 
	call readingobs(without_sat='yes')  ! спутников, то есть Ntype известно, но список самих спутников в заголовке не дан,  
else                                    ! то процедура readingobs также будет работать по другому пути. Поэтому и if-elseif-else 
	call readingobs(no_Nsat='yes')      !
endif           						!

write(*,*) 'Считывание данных rinex-файла завершено'

! Можно вывести интересующие данные на экран или в файл. Данные следует брать из этих массивов:
!
! list_time(1:Nobs) -- содержит дату и время каждого наблюдения (yy mm dd hh mm sec.00)
! datavalues(1:Ntype,1:Nsat,1:Nobs) -- содержит наблюдательные данные для каждого спутника для 
! 		каждого наблюдения (value LLI power, value2 LLI2 power2, ...)*количество спутников в
!		конкретном моменте наблюдений
! what_the_sats(1:Nsat,1:Nobs) -- содержит имена всех спутников в отдельно взятом наблюдении 
! all_numsat(1:Nobs) -- содержит количество имеющихся спутников с данными в отдельно взятом наблюдении

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!==================================================================================================
!**************************************************************************************************
!==================================================================================================
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

! Программа закончилась. Дальше идет попытка проверки правильности считывания данных.
! Большой rinex файл (200 Мб, ежесекундные наблюдения в течение 12 часов, 14 типов наблюдений,
! 55 спутников всего, в среднем по 17 в одном моменте наблюдений) был считан. И дальше написан код,
! заставляющий программу вывести в файл количество имеющихся наблюдательных данных каждого типа
! для каждого спутника. В том формате, в котором это дано в заголовке ( PRN / # OF OBS ).
! И, о счастье, после долгих мучений -- попыток вывести результаты нужным образом,
! программа не выдавала правильные результаты сразу по всем спутникам, файл с результатами наконец
! оказался идентичен тому, что написано в заголовке. Ура!
! Я молодец (если бы не так много времени на это убил, времени, которое я должен был тратить на другие дела)

allocate(megaprov(1:Ntype,1:Nsat),justsats(0:Nsat,1:Nobs))
allocate(datavalues0(1:Ntype,0:Nsat,1:Nobs))

!do k=1,Nobs
!	forall (i=all_numsat(k)+1:Nsat) datavalues(:,i,k)%value=0.0_pr
!end do

!datavalues0(:,1:Nsat,:)=datavalues

!forall (i=1:Nsat) justsats(i,:)=i
!!justsats=0
!
!do k=1,Nobs
!	do i=1,all_numsat(k)
!		do j=1,Nsat
!			if (what_the_sats(i,k)==sat_names(j)) then
!				if (i/=j) forall (s=i:j:j-i) justsats(s,k)=j+i-s
!			end if
!		end do
!	end do
!end do
!
!!write(*,*) 'forall ended'
!
!forall (k=1:Nobs) datavalues0(:,(/justsats(:,k)/),k)=datavalues0(:,:,k)

justsats=0

do k=1,Nobs
	do i=1,all_numsat(k)
		do j=1,Nsat
			if (what_the_sats(i,k)==sat_names(j)) then
				justsats(i,k)=j
			end if
		end do
	end do
end do

do k=1,Nobs
	forall (j=1:all_numsat(k),justsats(j,k)/=0) datavalues0(:,justsats(j,k),k)=datavalues(:,j,k)
end do

forall (i=1:Ntype, k=1:Nsat) 
	megaprov(i,k)=count(abs(datavalues0(i,k,:)%value)>eps)
end forall

open(100,file='prnobs.dat')
	do k=1,Nsat
		write(100,'(3x,a3,9i6)') sat_names(k), megaprov(1:9,k)
		if (Ntype>9) then
			do i=10,Ntype,9
				write(100,'(3x,3x,9i6)') megaprov(i:min(Ntype,i+8),k)
			end do
		end if
	end do
close(100)




end program getall