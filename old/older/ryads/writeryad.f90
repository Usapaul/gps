program ryad
! *** Программа выдает первые несколько членов степенного ряда ***

implicit none
integer, parameter :: mp=8
integer :: i, j, n
integer(mp) :: ch, zn, chfact, znfact ! Пояснения к переменным будут даны по ходу
integer(mp), allocatable :: chis(:,:), znam(:,:) ! Пояснения в коде
integer, dimension(1:2) :: stepen ! Сюда запишется 2 числа : a и b, когда переменная выглядит так: x**(a*k+b)
integer :: ncifr ! Эта переменная понадобится для расчета количества цифр знаменателя (красивый вывод результата)
character(5) :: testing
! Далее идут строковые переменные, которые нужны для их комбинации в format
character(7) :: f0='",i9'
character(30) :: f1='("k=",i2,", k-тый член:     '
character(7) :: f2=',"x^",i'
character(13) :: f3=',"/",i'
character(50) :: f4
character(5) :: fs='1'

!---------------------------------------
call getarg(1,testing)
if (testing=='test') then ! Если нужно быстро показать, что программа работает, то сюда можно заранее вписать данные
    allocate(chis(0:2,1),znam(0:2,1))      ! Пояснения ко всем массивам будут ниже
    stepen=(/2,1/); n=6
    chis(:,1)=(/2,1,-2/); znam(:,1)=(/3,1,0/)
    go to 88
endif
!---------------------------------------
write(*,*) '-----------------------'
write(*,*) 'Будут выданы значения коэффициентов и степеней'
write(*,*) 'Нескольких первых членов степенного ряда. Суммирование по k'
write(*,*) 'Введите общий член ряда, задавая в нем множители:'
write(*,*) '-----------------------'
write(*,*) 'Введите степень при x. Выглядит это так: x**(a*k+b)'
write(*,*) 'Введите эти числа a и b'
read(*,*) stepen(1), stepen(2)
!---------------------------------------
write(*,*) 'В различных произведениях участвуют множители вида:'
write(*,*) '(a*k+b), b**(a*k), (a*k+b)!. Найдите для каждого множителя числа a и b'
write(*,*) 'Если не удается, то преобразуйте множители: a**(b+c)=(a**b)*(a**c)'
write(*,*) '-----------------------'
!---------------------------------------
write(*,*) 'Сначала введите все множители в числителе. Cколько их?'
read(*,*) n; allocate(chis(0:2,n)); call givefun(chis) ! Массив содержит информацию о каждом множителе в числителе
write(*,*) 'Теперь введите все множители в знаменателе. Cколько их?'
read(*,*) n; allocate(znam(0:2,n)); call givefun(znam) ! -----"------ в знаменателе
!---------------------------------------
write(*,*) 'И, наконец, введите количество интересуемых первых членов "N"'
write(*,*) 'Будут выданы N+1 членов (от нуля до N):'
read(*,*) n

88 open(777,file='ress.dat') ! Файл с данными для программы, которая считает значение суммы первых членов при данном x
    write(777,*) n
do i=0,n
    ch=1; zn=1 ! Конечные значения числителя и знаменателя записывается в эти переменные (по умолчанию это 1)
    call rachetfun(i,chis,ch); call rachetfun(i,znam,zn); call nod(ch,zn)
    write(777,*) ch, stepen(1)*i+stepen(2), zn
    if (stepen(1)*i+stepen(2) >= 10 ) fs='2' ! При достижении степенью переменной значения 10, кол-во цифр увеличилось
select case(zn) ! Выбирается наиболее подходящий вывод, чтобы "все было красиво"
case(0); write(*,*) 'Деление на ноль. Проверьте правильность введенных данных' ! Если вдруг знаменатель нулевой
case(1,-1); if (zn==-1) ch=-ch
    select case(ch)
    case(0); write(*,f1//'         ",i1)') i, 0
    case(1); write(*,f1//'         "'//f2//fs//')') i, stepen(1)*i+stepen(2)
    case(-1); write(*,f1//'        -"'//f2//fs//')') i, stepen(1)*i+stepen(2)
    case default; write(*,f1//f0//f2//fs//')') i, ch, stepen(1)*i+stepen(2)
    end select
case default; if (zn<0) then; ch=-ch; zn=-zn; endif
    call kolvocifr(zn,ncifr); write(f4,*) ncifr
    select case(ch)
    case(0); write(*,f1//'         ",i1)') i, 0
    case(1); write(*,f1//'         "'//f2//fs//f3//f4//')') i, stepen(1)*i+stepen(2), zn
    case(-1); write(*,f1//'        -"'//f2//fs//f3//f4//')') i, stepen(1)*i+stepen(2), zn
    case default; write(*,f1//f0//f2//fs//f3//f4//')') i, ch, stepen(1)*i+stepen(2), zn
    end select
end select
enddo
   close(777)

contains

subroutine givefun(mass)
! *** Процедура получает от пользователя данные об общем члене ряда
implicit none
integer, parameter :: mp=8
integer(mp), dimension(0:,1:) :: mass
integer :: i
!---------------------------------------
do i=1,size(mass,dim=2)
    write(*,*) 'Если следующий множитель вида a*k+b, введите "1"'
    write(*,*) 'Если следующий множитель вида b^(a*k), введите "2"'
    write(*,*) 'Если следующий множитель вида (a*k+b)!, введите "3"'
    read(*,*) mass(0,i) ! Для каждого множителя записывается информация о его виде
    do while (mass(0,i)<1 .or. mass(0,i)>3) ! Программа потребует ввести правильное число
        write(*,*) 'Введено число, отличное от {1,2,3}. Введите заново'
        read(*,*) mass(0,i)
    enddo

    select case(mass(0,i)) ! Для каждого множителя теперь записываются числа a и b
    case(1); write(*,*) 'Множитель (a*k+b). Введите a и b'
        read(*,*) mass(1,i), mass(2,i)
    case(2); write(*,*) 'Множитель b^(a*k). Введите a и b'
        read(*,*) mass(1,i), mass(2,i)
    case(3); write(*,*) 'Множитель (a*k+b)!, введите a и b'
        read(*,*) mass(1,i), mass(2,i)
    end select
enddo

end subroutine givefun

subroutine rachetfun(i,mass,res)
! *** Процедура рассчитывает значения числителя и знаменателя, используя известные данные
implicit none
integer, parameter :: mp=8
integer, intent(in) :: i
integer(mp), dimension(0:,1:), intent(in) :: mass
integer(mp), intent(out) :: res
integer(mp) :: resfact
integer :: j

res=1
!---------------------------------------
do j=1,size(mass,dim=2) ! Умножаются все множители
    select case(mass(0,j))
    case(1); res=res*(mass(1,j)*i+mass(2,j))
    case(2); res=res*(mass(2,j)**(mass(1,j)*i))
    case(3); call factorial(mass(1:2,j),i,resfact); res=res*resfact ! Факториал считается в отдельной процедуре
    end select
enddo

end subroutine rachetfun

subroutine factorial(AandB,k,fact)
! *** Процедура считает факториал числа (a*k+b)
! *** Где AandB(1)=a, AandB(2)=b
implicit none
integer, parameter :: mp=8
integer(mp), dimension(1:2), intent(in) :: AandB ! Массив со значениями a и b для (a*k+b)!
integer, intent(in) :: k ! Число k для (a*k+b)
integer(mp), intent(out) :: fact ! Результат - факториал данного числа
integer :: i, a

a=AandB(1)*k+AandB(2); fact=1
!---------------------------------------
if (a>=2) then
do i=2,a
fact=fact*i
enddo; endif

end subroutine factorial

subroutine kolvocifr(a,n)
! *** Процедура выполняет подсчет количества цифр в данном числе
implicit none
integer, parameter :: mp=8
integer(mp), intent(in) :: a
integer(mp) :: b
integer, intent(out) :: n

b=a; n=0
!---------------------------------------
if (b<0) then
b=abs(b); n=1; endif
do while (b>0)
    b=(b-mod(b,10))/10; n=n+1
enddo

end subroutine kolvocifr

subroutine nod(ch,zn)
! *** Процедура выполняет сокращение неправильной дроби ch/zn
implicit none
integer, parameter :: mp=8
integer(mp) :: ch, zn, a, b, ost

a=max(abs(ch),abs(zn)); b=min(abs(ch),abs(zn))
ost=1
!---------------------------------------
if (b/=0) then
    do while (ost>0)
        ost=mod(a,b); a=b; b=ost
    enddo
    ch=ch/a; zn=zn/a
endif

end subroutine nod

end program ryad
