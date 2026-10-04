program fdtd2d
    implicit none
    ! Фундаментальные постоянные и параметры моделирования
    real(8), parameter::pi = acos(-1.0)
    real(8), parameter::cc = 2.99792458e8              !скорость света в свободном пространстве
    real(8), parameter::muz = 4.0*pi*1.0e-7            !магнитная постоянная
    real(8), parameter::epsz = 1.0/(cc*cc*muz)         !электрическая постоянная
    real(8), parameter::etaz=sqrt(muz/epsz)
    real(8), parameter::freq = 5.0e+9                  !центральная частота источника
    real(8), parameter::lambda = cc/freq               !центральная длина волны источника
    real(8), parameter::omega = 2.0*pi*freq
    ! Параметры сетки
    integer, parameter::ie = 300           !число ячеек сетки по оси x
    integer, parameter::je = 300            !число ячеек сетки по оси y
    integer, parameter::ib=ie+1
    integer, parameter::jb=je+1

    integer, parameter::is=150            !расположение источника в направлении оси z
    integer, parameter::js=150          !расположение источника в направлении оси z
    ! Числовые неконстантные переменные
    !dx-пространственный шаг квадратной решетки, dt-временной шаг
    real(8) dx, dt
    !параметры среды
    real(8) eps, mur
    !коэффициенты обновления
    real(8) c, d
    !nt-общее число временных шагов, n-итератор
    integer nt, n
    ! Массивы для полей
    real(8) ex(ie,jb)           !поля на основной решетке
    real(8) ey(ib,je)
    real(8) hz(ie,je)

    ! Динамический массив для источника
    real(8), allocatable::source(:)
    ! Предыдущие значения полей для граничных условий Мура — поля на границах с предыдущего шага
    real(8) exprev1(ie), exprev2(ie), exprevjb(ie), exprevjbm(ie)
    real(8) eyprev1(je), eyprev2(je), eyprevib(je), eyprevibm(je) 
    ! Коэффициент для граничных условий
    real(8) rbc
    ! Массивы для сохранения результатов
    real(8) exres(ie, jb), eyres(ib, je), hzres(ie, je)
    real(8) tres
    ! Метка для файла сохранения и временные отсчеты
    integer results, countstart, countend, countrate, countmax 
    real(8) elapsed
    character(len=32) :: arg

    ! Массив для замены цикла инициализации источника
    real(8), allocatable::t(:)

    print *, "Source wavelength=", lambda

    ex = 0.0; ey = 0.0; hz = 0.0
    dx=3.0e-3
    dt=dx/(2.0*cc)  
    rbc=(cc*dt - dx)/(cc*dt + dx)
    nt=300 
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) nt
    end if
    eps=1.0
    mur=1.0

    ! Излучение волны
    allocate(source(nt), t(nt))
    source=0.0  
    t = [(n*dt, n = 1, nt)]
    source=sin(omega*t)

    ! Коэффициенты обновления
    c=dt/epsz/dx
    d=dt/muz/dx

    ! НАЧАЛО ЦИКЛА РАСЧЕТА РАСПРОСТРАНЕНИЯ
    call system_clock(countstart, countrate, countmax)
    do n=1, nt
        ! Сохранение граничных значений электрических полей перед обновлением
        exprev1=ex(:, 1); exprev2=ex(:, 2); exprevjb=ex(:, jb); exprevjbm=ex(:, jb-1)
        eyprev1=ey(1, :); eyprev2=ey(2, :); eyprevib=ey(ib, :); eyprevibm=ey(ib-1, :)

        ! Обновление электрического поля (EX и EY)
        ex(:,2:je)=ex(:,2:je)+c*(hz(:,2:je)-hz(:,1:je-1))
        ey(2:ie,:)=ey(2:ie,:)+c*(hz(1:ie-1,:)-hz(2:ie,:))

        ! Граничные условия Мура для EX
        ex(:, 1)=exprev2+rbc*(ex(:, 2)-exprev1)
        ex(:, jb)=exprevjbm+rbc*(ex(:, jb-1)-exprevjb)

        ! Граничные условия Мура для EY 
        ey(1, :)=eyprev2+rbc*(ey(2, :)-eyprev1)
        ey(ib, :)=eyprevibm+rbc*(ey(ib-1, :)-eyprevib)

        ! Обновление магнитного поля (HZ)
        hz(1:ie,1:je)=hz(1:ie,1:je)+d*(ex(1:ie,2:jb)-ex(1:ie,1:je)+ &
        & ey(1:ie,1:je)-ey(2:ib,1:je))
        hz(is,js)=source(n)
    end do
    call system_clock(countend, countrate, countmax)
    elapsed = real(countend - countstart) / real(countrate)
    print *, 'Calculation time = ', elapsed, ' s'

    ! Сохранение полей для визуализации
    exres = ex
    eyres = ey
    hzres = hz
    tres = dt * nt / 1.0e-9

    open(newunit=results, file='results/fields2d.bin', &
    form='unformatted', access='stream', status='replace')
    write(results) int(ie, 8), int(jb, 8), int(ib, 8), int(je, 8), int(1, 8)
    write(results) real(dx, 8), real(dt, 8)
    write(results) real(tres, 8)
    write(results) real(exres, 8)
    write(results) real(eyres, 8)
    write(results) real(hzres, 8)
    close(results)

end program fdtd2d