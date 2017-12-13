program time
implicit none

real(8), allocatable :: ryad(:), poss(:,:), tx(:)
integer(8) :: i, k=1e6, n, j, m, kn
integer(8), allocatable :: pss(:), t(:), sort(:), sort1(:), average(:)
!-----------------------------
open(100,file='possbl.dat')
    read(100,'(2x,i8)') n
    allocate(poss(2,n),pss(n))
    do i=1,n
        read(100,*) poss(1:2,i)
    enddo
    pss(1:n)=idint(poss(2,1:n)*k/n)
close(100)
!-----------------------------
kn=sum(pss)
allocate(ryad(kn),t(kn),tx(kn))
m=0; j=0
do i=1,n
    j=pss(i)
    ryad(m+1:m+j)=poss(1,i)
    m=m+j
enddo
!-----------------------------
call random_number(tx)
t=idint(tx*kn)
allocate(sort(kn),sort1(kn),average(kn-1))
forall (i=1:kn) sort1(i)=ryad(t(i)+1)*3262+300*i

write(*,*) '    number of SN: ', kn

do i=1,kn
    sort(i)=minval(sort1)
    sort1(minloc(sort1))=maxval(sort1)+1
    if (mod(i,kn/10)==0) write(*,*) i*100/kn, '%'
enddo
!-----------------------------
open(777,file='time.dat')
    do i=1,kn-1
        write(777,*) sort(i+1)-sort(i)
        average(i)=sort(i+1)-sort(i)
    enddo
write(*,*) '     average: ', sum(average)/kn
close(777)


end program time