program fdtd2d
    implicit none
    ! Фундаментальные постоянные и параметры моделирования
    real(8), parameter :: pi = acos(-1.0)
    real(8), parameter :: cc = 2.99792458e8                ! скорость света в свободном пространстве
    real(8), parameter :: muz = 4.0*pi*1.0e-7              ! магнитная постоянная
    real(8), parameter :: epsz = 1.0/(cc*cc*muz)           ! электрическая постоянная
    real(8), parameter :: etaz = sqrt(muz/epsz)
    real(8), parameter :: freq = 5.0e+9                   ! центральная частота источника
    real(8), parameter :: lambda = cc/freq                ! центральная длина волны источника
    real(8), parameter :: omega = 2.0*pi*freq

    ! Параметры сетки
    integer, parameter :: ie = 300           ! число ячеек сетки по оси x
    integer, parameter :: je = 300           ! число ячеек сетки по оси y
    integer, parameter :: ib = ie + 1
    integer, parameter :: jb = je + 1

    integer, parameter :: is = 150           ! расположение источника в направлении оси z
    integer, parameter :: js = 150           ! расположение источника в направлении оси z

    ! Числовые неконстантные переменные
    real(8) dx, dt                          ! dx — пространственный шаг, dt — временной шаг
    real(8) eps, mur                        ! параметры среды
    real(8) c, d                            ! коэффициенты обновления
    integer nt, n, i, j                     ! nt — общее число временных шагов, n, i, j — итераторы

    ! Массивы для полей
    real(8) ex(ie, jb)                      ! поля на основной решётке
    real(8) ey(ib, je)
    real(8) hz(ie, je)

    ! Динамический массив для источника
    real(8), allocatable :: source(:)
    real(8), allocatable :: t(:)

    ! Предыдущие значения полей для граничных условий Мура
    real(8) exprev1(ie), exprev2(ie), exprevjb(ie), exprevjbm(ie)
    real(8) eyprev1(je), eyprev2(je), eyprevib(je), eyprevibm(je)

    real(8) rbc                             ! коэффициент для граничных условий

    ! Массивы для сохранения результатов
    real(8) exres(ie, jb), eyres(ib, je), hzres(ie, je)
    real(8) tres
    ! Метка для файла сохранения и временные отсчёты
    integer results, countstart, countend, countrate, countmax
    real(8) elapsed
    character(len=32) :: arg

    print *, "Source wavelength = ", lambda

    ex = 0.0; ey = 0.0; hz = 0.0
    
    dx = 3.0e-3
    dt = dx / (2.0 * cc)
    rbc = (cc*dt - dx) / (cc*dt + dx)
    nt=300 
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) nt
    end if
    eps = 1.0
    mur = 1.0

    allocate(source(nt), t(nt))
    source = 0.0
    do n = 1, nt
        source(n) = sin(omega*n*dt)
    end do
    c = dt / epsz / dx
    d = dt / muz / dx

    call system_clock(countstart, countrate, countmax)

    do n = 1, nt
        ! Сохранение граничных значений электрических полей перед обновлением
        ! EX
        do i = 1, ie
            exprev1(i) = ex(i, 1)
            exprev2(i) = ex(i, 2)
            exprevjb(i) = ex(i, jb)
            exprevjbm(i) = ex(i, jb-1)
        end do

        ! EY
        do j = 1, je
            eyprev1(j) = ey(1, j)
            eyprev2(j) = ey(2, j)
            eyprevib(j) = ey(ib, j)
            eyprevibm(j) = ey(ib-1, j)
        end do

        ! Обновление электрического поля EX 
        do i = 1, ie
            do j = 2, je
                ex(i, j) = ex(i, j) + c * (hz(i, j) - hz(i, j-1))
            end do
        end do

        ! Обновление электрического поля EY
        do i = 2, ie
            do j = 1, je
                ey(i, j) = ey(i, j) + c * (hz(i-1, j) - hz(i, j))
            end do
        end do

        ! Граничные условия Мура для EX
        do i = 1, ie
            ex(i, 1) = exprev2(i) + rbc * (ex(i, 2) - exprev1(i))
            ex(i, jb) = exprevjbm(i) + rbc * (ex(i, jb-1) - exprevjb(i))
        end do

        ! Граничные условия Мура для EY
        do j = 1, je
            ey(1, j) = eyprev2(j) + rbc * (ey(2, j) - eyprev1(j))
            ey(ib, j) = eyprevibm(j) + rbc * (ey(ib-1, j) - eyprevib(j))
        end do

        ! Обновление магнитного поля HZ
        do i = 1, ie
            do j = 1, je
                hz(i, j) = hz(i, j) + d * (ex(i, j+1) - ex(i, j) + ey(i, j) - ey(i+1, j))
            end do
        end do

        ! Источник
        hz(is, js) = source(n)
    end do

    call system_clock(countend, countrate, countmax)
    elapsed = real(countend - countstart) / real(countrate)
    print *, 'Calculation time = ', elapsed, ' s'
    ! Сохранение полей для визуализации
    exres = ex
    eyres = ey
    hzres = hz
    tres = dt * nt / 1.0e-9

    open(newunit=results, file='results/fields2d_scalar.bin', &
    form='unformatted', access='stream', status='replace')
    write(results) int(ie, 8), int(jb, 8), int(ib, 8), int(je, 8), int(1, 8)
    write(results) real(dx, 8), real(dt, 8)
    write(results) real(tres, 8)
    write(results) real(exres, 8)
    write(results) real(eyres, 8)
    write(results) real(hzres, 8)
    close(results)

end program fdtd2d
