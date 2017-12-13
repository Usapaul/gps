program summa
! *** Программа рассчитывает сумму заданного ряда при данном x

implicit none
integer, parameter :: mp=8
real(mp) :: res, x ! Искомая сумма - res
integer(mp) :: ch, zn ! Числитель и знаменатель
integer :: stepen, i, n ! Степень икса - stepen
character(18) :: fx, fxp ! Переменные для оператора format

res=0
!---------------------------------------
open(100,file='ress.dat')
    read(100,*) n ! Читается количество интересующих членов
    write(*,'("Введите x, для которого считается сумма членов ряда от 0 до ",i2)') n
    read(*,*) x ! Для этого x и будет вычислена сумма конечного ряда
    do i=0,n
	read(100,*) ch, stepen, zn
	res=res+ch*x**(stepen)/zn
    enddo
close(100)
!---------------------------------------
write(fx,*) mp*2+1  ! Количество цифр, отведенных на x
write(fxp,*) mp*2-1 ! Количество цифр, отведенных на дробную часть x
write(*,'("Сумма членов ряда от 0 до ",i2," при x=",f'//fx//'.'//fxp//')') n, x
write(*,*) res

end program summa
