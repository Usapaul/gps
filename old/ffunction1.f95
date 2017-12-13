module ffunction1

contains

function Ffun1(X) result(Ff)
! *** Функция Ffun(X) получает вектор X, и возвращает вектор Ffun той же размерности
implicit none
real(8), dimension(1:), intent(in) :: X
real(8), dimension(1:size(X)) :: Ff
integer :: i, n
real(8), dimension(1:size(X),1:size(X)) :: XYZP ! 1 спутник - 1 строка

n=size(X)
!---------------------------------------
open(500,file='XYZPsat.dat')
    read(500,*) XYZP(:,1)
    read(500,*) XYZP(:,2)
    read(500,*) XYZP(:,3)
    read(500,*) XYZP(:,4)
close(500)

Ff=(XYZP(1,:)-X(1))**2+(XYZP(2,:)-X(2))**2+(XYZP(3,:)-X(3))**2-(XYZP(4,:)+X(4))**2

end function Ffun1

end module ffunction1