program graphSN
implicit none

real(8) :: r, RG, R0
real(8) :: h, fi
real(8), parameter :: pi=3.14159265359
integer :: n=10000
!-----------------------------
RG=15.0d0 !кпк
R0=8.0d0 !кпк
h=(RG+R0)/n
!-----------------------------
open(777,file='integral.dat')
write(777,'("# ",i8)') n
r=0.0d0
write(777,*) r, 0.0d0
do while (r+h<RG-R0)
    r=r+h
    write(777,*) r, (r/RG)**2
enddo
!-----------------------------
do while (r+h<sqrt(RG**2-R0**2))
    r=r+h
    fi=acos((RG**2+R0**2-r**2)/2/RG/R0)
    write(777,*) r, (r/RG)**2*(1.0d0-asin(sin(fi)*RG/r)/pi)+fi/pi-R0/RG*sin(fi)/pi
enddo
!-----------------------------
do while (r+h<RG+R0)
    r=r+h
    fi=acos((RG**2+R0**2-r**2)/2/RG/R0)
    write(777,*) r, (r/RG)**2*asin(sin(fi)*RG/r)/pi+fi/pi-R0/RG*sin(fi)/pi
enddo
!-----------------------------
write(777,*) RG+R0, 1.0d0
close(777)

end program graphSN

