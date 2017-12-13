module readingheadermod
use typeconst

implicit none

character(80), dimension(:), allocatable :: header_array
character(2),  dimension(:), allocatable :: type_names
character(3),  dimension(:), allocatable :: sat_names

contains

subroutine readingheader()
	implicit none

	integer, parameter :: maxNheaderlines=5000, maxNtype=50
	character(80), dimension(:), allocatable :: header_array_dump
	character(2),  dimension(:), allocatable :: type_names_dump
	integer :: i, n=0, k, j

	Ntype=0; Nsat=0; Nheaderlines=0
	allocate(header_array_dump(0:maxNheaderlines),type_names_dump(maxNtype))
	!----------------------------------
	open(100,file='rinexfortran.dat',status='old')
		i=0
		header_array_dump(i)=' '
		do while (header_array_dump(i)(61:80)/='END OF HEADER       '.and.i<size(header_array_dump))
			! Считываются отдельно метки заголовка и данные в них
			i=i+1
			read(100,'(a80)') header_array_dump(i)
		end do
		if (i==size(header_array_dump)) then
			stop 'END OF HEADER не найден!'
		end if

		n=i ! Здесь n равно количеству строк в заголовке
		Nheaderlines=n ! Главная программа тоже должна знать, сколько в заголовке строк	
	close(100,status='keep')
	allocate(header_array(n))              ! Этот массив не занимает лишнее место,
	header_array(1:n)=header_array_dump(1:n)
	deallocate(header_array_dump)            ! по сравнению с этим, который называется dump

	k=1
	do while (k<=n)
    select case(header_array(k)(61:80))
    	case('APPROX POSITION XYZ ')
    		read(header_array(k),'(3(f14.4))') approx_XYZ(1:3)
    		k=k+1
		case('# / TYPES OF OBSERV ')
			read(header_array(k),'(i6,*(4x,a2))') Ntype, type_names_dump(1:min(Ntype,9))
			if (Ntype<=0) stop 'Количество типов наблюдений указано некорректно'
			! Но типов может быть больше 9, тогда запись продолжается в следующих строках (по 9 шт)
			if (Ntype>9.and.Ntype<=maxNtype) then
				do i=10,Ntype,9
					k=k+1
					read(header_array(k),'(6x,*(4x,a2))') type_names_dump(i:min(Ntype,i+8))
				end do
			else if (Ntype>=maxNtype) then  ! Я тут указал maxNtype=50, должно хватать, так что это просто формальность.
				stop 'Ntype > max N of types in the program'
			end if
			allocate(type_names(Ntype))
			type_names(1:Ntype)=type_names_dump(1:Ntype)
			deallocate(type_names_dump)
			k=k+1
		case('# OF SATELLITES     ')
    		read(header_array(k),'(i6)') Nsat
    		k=k+1
    	case('PRN / # OF OBS      ')
    		if (Nsat>0) then
    			allocate (sat_names(Nsat))
    			k=k-1
    			do i=1,Nsat
    				k=k+1
    				read(header_array(k),'(3x,a3)') sat_names(i)
    				if (Ntype>9) then
    					do j=2,(Ntype-1)/9+1
    						k=k+1
    					end do
    				end if
    			end do
    		else
    			continue
    		end if
    		k=k+1
		case('RINEX VERSION / TYPE','PGM / RUN BY / DATE ','COMMENT             ')
    		k=k+1
    	case('MARKER NAME         ','MARKER NUMBER       ','OBSERVER / AGENCY   ')
    		k=k+1
		case('REC # / TYPE / VERS ','ANT # / TYPE        ','ANTENNA: DELTA H/E/N')
			k=k+1
    	case('WAVELENGTH FACT L1/2')
    		k=k+1
    	case('INTERVAL            ','TIME OF FIRST OBS   ','TIME OF LAST OBS    ')
    		k=k+1
    	case('RCV CLOCK OFFS APPL ','LEAP SECONDS        ')
			k=k+1    	
    	case('END OF HEADER       ')
    		k=k+1
    	case default
    		write(*,*) ' Обнаружена неопознанная метка в заголовке:', header_array(k)(61:80)
    		k=k+1
    end select
	end do


end subroutine readingheader

end module readingheadermod