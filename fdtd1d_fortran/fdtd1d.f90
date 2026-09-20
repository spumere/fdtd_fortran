program fdtd1d
    implicit none
    ! Фундаментальные постоянные и параметры материала
    real, parameter::pi = acos(-1.0)
    real, parameter::cc = 2.99792458e8              !скорость света в свободном пространстве
    real, parameter::muz = 4.0*pi*1.0e-7            !магнитная постоянная
    real, parameter::epsz = 1.0/(cc*cc*muz)         !электрическая постоянная
    real, parameter::freq = 1.0e+9                  !частота источника
    real, parameter::lambda = cc/freq               !длина волны источника
    real, parameter::omega = 2.0*pi*freq
    real, parameter::eps=1.0, sig=5.0e-3
    ! Параметры сетки, массивы полей и 
    ! коэффициенты для пространственной области с непроницаемой средой
    integer, parameter::ie = 200, ih = ie-1         !число отсчетов Hy и Ez на сетке
    real dx, dt, omegadt, x(ih), ez(ie), hy(ih), scfact, ca, cb, rbc, elapsed, alpha, beta, tgdelta, eta_ez, eta_hy, xi(ih)
    complex zc !характеристическое сопротивление для аналитического решения
    integer::i, n, nmax, results, count_start, count_end, count_rate, count_max 
    !динамические массивы для хранения результатов численного и аналитического решений      
    real, allocatable :: ez_res(:,:), hy_res(:,:), t_res(:), ez_an(:,:), hy_an(:,:) 
    dx = lambda/20.0                                !пространственный шаг одномерной решетки
    dt = dx/cc                                      !временной шаг
    omegadt = omega*dt
    nmax = nint(12.0e-9/dt)                         !общее число временных шагов
    x = [(i*dx, i = 1, ih)]
    scfact=dt/muz/dx;
    ca=(1.0-(dt*sig)/(2.0*epsz*eps))/(1.0+(dt*sig)/(2.0*epsz*eps));
    cb=scfact*(dt/epsz/eps/dx)/(1.0+(dt*sig)/(2.0*epsz*eps));
    ez = 0.0
    hy = 0.0
    allocate(ez_res(ie, nmax), hy_res(ih, nmax), t_res(nmax))
    ! НАЧАЛО ЦИКЛА РАСЧЕТА РАСПРОСТРАНЕНИЯ
    call system_clock(count_start, count_rate, count_max)
    do n = 1, nmax
        ez(1)=scfact*sin(omegadt*n)
        rbc=ez(ih)
        do i = 2, ih
            ez(i)=ca*ez(i)+cb*(hy(i)-hy(i-1))
        end do
        ez(ie)=rbc
        do i = 1, ih
            hy(i)=hy(i)+ez(i+1)-ez(i)
        end do
        ! Запись полей в массивы
        do i = 1, ie
            ez_res(i, n)=ez(i)/scfact
        end do
        do i = 1, ih
            hy_res(i, n)=hy(i)
        end do
        t_res(n)=dt*n/1.0e-9
    end do
    call system_clock(count_end, count_rate, count_max)

    ! Вычисление аналитического решения по Пименов и др. Техническая электродинамика. М.: Радио и связь, 2000.
    allocate(ez_an(ie, nmax), hy_an(ih, nmax))
    tgdelta = sig/omega/epsz
    alpha = omega*sqrt(epsz*muz/2*(sqrt(1+tgdelta**2)-1))
    beta = omega*sqrt(epsz*muz/2*(sqrt(1+tgdelta**2)+1))
    zc = sqrt(muz/epsz/(1.0-(0.0, 1.0)*tgdelta))

    xi = x - dx
    
    do n = 1, nmax
        do i = 1, ih
            if (real(n)*dt >= xi(i)/cc) then
                !волна принимает ненулевое значение только в тех точках, до которых физически успела дойти
                ez_an(i, n) = exp(-alpha*xi(i)) * &
                            sin(omegadt*n - beta*xi(i))
                hy_an(i, n) = -1.0/abs(zc) * exp(-alpha*x(i)) * &
                            sin(omegadt*n - beta*xi(i) - atan(tgdelta)/2.0)
            else
                ez_an(i, n) = 0.0
                hy_an(i, n) = 0.0
            end if
        end do
    end do

    elapsed = real(count_end - count_start) / real(count_rate)
    print *, 'Calculation time = ', elapsed, ' s'

    ! Вычисление интеграла перекрытия (нормированнного скалярного произведения) аналитического и численного решений
    ! На всем временном промежутке распространения
    eta_ez = sum(ez_res*ez_an)/sqrt(sum(ez_res*ez_res))/sqrt(sum(ez_an*ez_an))
    eta_hy = sum(hy_res*hy_an)/sqrt(sum(hy_res*hy_res))/sqrt(sum(hy_an*hy_an))

    print *, 'Ez overlap across entire task = ', eta_ez
    print *, 'Hy overlap across entire task = ', eta_hy

    ! Сохранение массивов в файл и построение графиков в Python
    open(newunit=results, file='results/fields.bin', &
        form='unformatted', access='stream', status='replace')
    write(results) int(ih, 4), int(nmax, 4)
    write(results) real(x(1:ih), 4)
    write(results) real(t_res(1:nmax), 4)
    do n = 1, nmax
        write(results) real(ez_res(1:ih, n), 4)
    end do
    do n = 1, nmax
        write(results) real(ez_an(1:ih, n), 4)
    end do
    do n = 1, nmax
        write(results) real(hy_res(1:ih, n), 4) 
    end do
    do n = 1, nmax
        write(results) real(hy_an(1:ih, n), 4)  
    end do
    close(results)
    deallocate(ez_res, hy_res, t_res, ez_an, hy_an)
    call execute_command_line('python results/plot.py')

end program fdtd1d