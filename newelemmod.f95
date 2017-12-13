module newelemmod

implicit none

!interface newelem
!	logical function newelem(array,new)
!		character(*), dimension(1:), intent(in) :: array
!		character(len(array(1))),    intent(in) :: new
!	end function newelem
!
!	newelem, differentelem
!end interface newelem

interface newelem
	module procedure newelem, differentelem
end interface

contains

logical function newelem(array,new)
	character(*), dimension(1:), intent(in) :: array
	character(len(array(1))),    intent(in) :: new
	newelem=.NOT.ANY(array==new)
end function newelem

logical function differentelem(old,new)
	character(*),        intent(in) :: old
	character(len(old)), intent(in) :: new
	differentelem=old/=new
end function differentelem
!--------------------------------------
!--------------------------------------
!--------------------------------------
subroutine readdump(idfile,nlines)
	! Процедура перемещает файловый указатель на nlines строк вперед
	implicit none

	integer, intent(in) :: nlines, idfile
	character(1) :: dump
	character(10) :: nlineschar

	write(nlineschar,'(i8)') nlines-1
	read(idfile,'('//nlineschar//'/,a1)') dump

end subroutine readdump

end module newelemmod