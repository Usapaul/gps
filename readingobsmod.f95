module readingobsmod
use typeconst
use readingheadermod
use newelemmod

implicit none

integer :: nlinesforsat ! В переменной указано, сколько строк отводится под запись набл. данных для одного спутника
type(datewithsec), dimension(:),     allocatable :: list_time     ! Здесь хранится дата и время для каждого наблюдения
type(obs_data),    dimension(:,:,:), allocatable :: datavalues    ! а здесь все значения наблюдательных данных
character(3),      dimension(:,:),   allocatable :: what_the_sats ! а вот здесь хранятся списки имен спутников,
integer,           dimension(:),     allocatable :: all_numsat    ! а тут -- сколько спутников было в каждом наблюдении

contains 

subroutine readingobs(without_sat,no_Nsat)
	implicit none

	integer, parameter :: maxNobs=100000
	character(3), intent(in), optional :: without_sat, no_Nsat
	type(obs_firstline) :: line1
	type(datewithsec), dimension(maxNobs)         :: list_time_dump 
	type(obs_data), dimension(:,:,:), allocatable :: datavalues_dump
	character(3),   dimension(:,:),   allocatable :: snames
	integer,        dimension(maxNobs)            :: all_numsat_dump
	integer :: i
	 									       ! Если количество спутников numsat<12, то надо сместить файловый
	character(3), dimension(12) :: snames_dump ! указатель на (12-numsat) позиций, записав сюда пробелы из
											   ! первой строки до самого delta_t 
	
	nlinesforsat=(Ntype-1)/5+1
	Nobs=0
	!--------------------------------------	
	if (present(no_Nsat)) then
		call wantsatnames('no_Nsat ')
	else if (present(without_sat)) then
		call wantsatnames('no_names')
	else
		continue ! Значит, количество спутников и их имена были получены из заголовка
	end if
	allocate(datavalues_dump(1:Ntype,1:Nsat,1:maxNobs),snames(1:Nsat,1:maxNobs))
	open(100,file='rinexfortran.dat',status='old')
		call readdump(idfile=100,nlines=Nheaderlines) ! Процедура перемещает файловый указатель на nlines строк
		wholefile: do
			read(100,1111,end=999) line1%datetime%datehm, line1%datetime%sec, line1%flag, line1%numsat,&
				& snames(1:min(line1%numsat,12),Nobs+1), snames_dump(min(line1%numsat,12)+1:12), line1%delta_t
			1111 format (1x,i2.2,4(1x,i2),f11.7,2x,i1,i3,12(a3),f12.9)
			if (line1%numsat>12) then
				do i=13,line1%numsat,12
					read(100,'(32x,12(a3))') snames(i:min(line1%numsat,i+11),Nobs+1)
				end do
			end if
			Nobs=Nobs+1
			list_time_dump(Nobs)=line1%datetime
			all_numsat_dump(Nobs)=line1%numsat
			call readoneobs(line1%numsat,datavalues_dump(1:Ntype,1:line1%numsat,Nobs))
			if (mod(Nobs,500)==0) write(*,*) ' Nobs= ', Nobs ! Чтобы видеть процесс, когда наблюдений дофига
		end do wholefile
		999 continue
		allocate(datavalues(Ntype,Nsat,Nobs),list_time(Nobs),what_the_sats(Nsat,Nobs),all_numsat(Nobs))
		datavalues   =datavalues_dump(:,:,1:Nobs)
		list_time    =list_time_dump(1:Nobs)
		what_the_sats=snames(:,1:Nobs)
		all_numsat   =all_numsat_dump(1:Nobs)
	close(100,status='delete')
end subroutine readingobs
!======================================
subroutine readoneobs(numsat,values)
	implicit none

	integer, intent(in) :: numsat
	type(obs_data), dimension(1:,1:), intent(inout) :: values
	integer :: i, j, k

	!----------------------------------
	do k=1,numsat
		do j=1,nlinesforsat
			read(100,'(5(f14.3,i1,i1))') (values(i,k)%value, values(i,k)%LLI, values(i,k)%power, i=(j-1)*5+1,min(Ntype,j*5))
		end do
	end do
end subroutine readoneobs
!======================================
subroutine wantsatnames(whatsapp)
	implicit none

	integer, parameter :: maxNsat=1000
	character(8), intent(in) :: whatsapp
	integer :: numsat=0, i=0, localNsat=0
	character(3), dimension(maxNsat) :: snames, allnames ! Имена в одном наблюдении и список уникальных имен
	!--------------------------------------
	open(100,file='rinexfortran.dat',status='old')
		call readdump(idfile=100,nlines=Nheaderlines)
		read(100,'(29x,i3,12(a3))') numsat, snames(1:min(numsat,12))
			if (numsat>12) then
				do i=13,numsat,12
					read(100,'(32x,12(a3))') snames(i:min(numsat,i+11))
				end do
			end if
			allnames(1:numsat)=snames(1:numsat)
			localNsat=numsat
			call readdump(idfile=100,nlines=numsat*nlinesforsat)
			
			! Если в первом наблюдении оказались все Nsat спутников, то процедура завершается через return:
			if (whatsapp=='no_names'.and.localNsat==Nsat) then
				allocate(sat_names(Nsat))
				sat_names(1:Nsat)=allnames(1:Nsat)
				close(100,status='keep')
				return
			end if
			!--------------------------
			! Первый прогон завершается тем, что были получены имена спутников в первом
			! наблюдении. Если не было известно Nsat, то счетчик уникальных имен начинается
			! с первого прогона, и дальше localNsat будет увеличиваться при встрече новых имен спутников
			! А когда на каком-то шаге localNsat достигнет Nsat, массив sat_names будет заполнен 
			! полученными к этому шагу именами спутников, хранящихся в массиве allnames(1:localNsat)
			! ---
			! Если число спутников Nsat не было известно после считывания заголовка, то придется его
			! получить в этой процедуре. Здесь производится считывание всего файла до конца с сохранением
			! в массив allnames(localNsat) на каждом шаге новых названий спутников
			! (и тогда localNsat увеличивается на единицу), а в конце Nsat присваевается localNsat
			!--------------------------
		everyobs: do
			read(100,'(29x,i3,12(a3))',end=999) numsat, snames(1:min(numsat,12))
			if (numsat>12) then
				do i=13,numsat,12
					read(100,'(32x,12(a3))') snames(i:min(numsat,i+11))
				end do
			end if
			! Запись нового встретившегося имени происходит вот так:
			do i=1,numsat
				if (newelem(allnames(1:localNsat),snames(i)).eqv..TRUE.) then
					localNsat=localNsat+1
					allnames(localNsat)=snames(i)
					if (localNsat==Nsat) exit everyobs
				end if
			end do
			call readdump(idfile=100,nlines=numsat*nlinesforsat)
		end do everyobs
		! При достижении конца файла цикл everyobs do просто завершается через метку 999,
		! но может произойти выход из цикла everyobs раньше, если было известно Nsat, и на
		! каком-то шаге localNsat достигло значения Nsat
		999 continue
		if (localNsat==Nsat.and.whatsapp=='no_names') then
			continue ! Значит, все прошло хорошо
		else if (localNsat/=Nsat.and.whatsapp=='no_names') then
			stop 'localNsat/=Nsat'
		else if (whatsapp=='no_Nsat ') then
			Nsat=localNsat
		end if
		allocate(sat_names(1:Nsat))
		sat_names(1:Nsat)=allnames(1:Nsat)
	close(100,status='keep')
end subroutine wantsatnames

end module readingobsmod