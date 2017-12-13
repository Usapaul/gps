module ffunction

contains

function Ffun(X) result(Ff)
! *** Функция Ffun(X) получает вектор X, и возвращает вектор Ffun той же размерности
implicit none
real(8), dimension(1:), intent(in) :: X
real(8), dimension(1:size(X)) :: Ff
integer :: i, n
real(8) :: e, Man

n=size(X)
!---------------------------------------
open(400,file='eMan.dat')
read(400,*) e
read(400,*) Man
close(400)

Ff=0
Ff=X-e*sin(X)-Man

end function Ffun

end module ffunction