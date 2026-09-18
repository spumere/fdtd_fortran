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
    real dx, dt, omegadt, x(ih), ez(ie), hy(ih), scfact, ca, cb, rbc, elapsed
    integer::i, n, nmax, u_ez, u_hy, count_start, count_end, count_rate, count_max         
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
    open(newunit=u_ez, file='results/ez.dat', status='replace')
    open(newunit=u_hy, file='results/hy.dat', status='replace')
    ! НАЧАЛО ЦИКЛА ВО ВРЕМЕНИ
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
        ! Запись полей в файл для визуализации
        write(u_ez, '(A,F10.4)') '# t = ', n*dt/1.0e-9
        write(u_hy, '(A,F10.4)') '# t = ', n*dt/1.0e-9
        do i = 1, ih
            write(u_ez, '(2ES15.6)') x(i), ez(i)/scfact
            write(u_hy, '(2ES15.6)') x(i), hy(i)
        end do
        write(u_ez, '(/)')
        write(u_hy, '(/)')
    end do
    call system_clock(count_end, count_rate, count_max)
    elapsed = real(count_end - count_start) / real(count_rate)
    print '(A,F10.4,A)', 'Wall time = ', elapsed, ' s'
    close(u_ez)
    close(u_hy)
    call execute_command_line('gnuplot results\plot.gp')

end program fdtd1d