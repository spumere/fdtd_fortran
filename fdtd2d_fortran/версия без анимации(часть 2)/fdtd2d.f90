program fdtd2d
    implicit none
    ! Фундаментальные постоянные и параметры моделирования
    real, parameter::pi = acos(-1.0)
    real, parameter::cc = 2.99792458e8              !скорость света в свободном пространстве
    real, parameter::muz = 4.0*pi*1.0e-7            !магнитная постоянная
    real, parameter::epsz = 1.0/(cc*cc*muz)         !электрическая постоянная
    real, parameter::etaz=sqrt(muz/epsz)
    real, parameter::freq = 5.0e+9                  !центральная частота источника
    real, parameter::lambda = cc/freq               !центральная длина волны источника
    real, parameter::omega = 2.0*pi*freq
    ! Параметры сетки
    integer, parameter::ie = 300           !число ячеек сетки по оси x
    integer, parameter::je = 300           !число ячеек сетки по оси y
    integer, parameter::ib=ie+1
    integer, parameter::jb=je+1

    integer, parameter::is=150           !расположение источника в направлении оси z
    integer, parameter::js=150          !расположение источника в направлении оси z
    integer, parameter::iebc=8           !толщина левой и правой областей PML
    integer, parameter::jebc=8           !толщина передней и задней областей PML
    integer, parameter::orderbc=2
    integer, parameter::ibbc=iebc+1
    integer, parameter::jbbc=jebc+1
    integer, parameter::iefbc=ie+2*iebc
    integer, parameter::jefbc=je+2*jebc
    integer, parameter::ibfbc=iefbc+1
    integer, parameter::jbfbc=jefbc+1
    ! Числовые неконстантные переменные
    !dx-пространственный шаг квадратной решетки, dt-временной шаг
    real dx, dt, rmax 
    !переменные для определения источника и коэффициентов обновления
    real rtau, tau, delay, eaf, haf
    !параметры материалов
    real eps(2), sig(2), mur(2), sim(2)
    !nmax-общее число временных шагов, media-число сред, i, j, n-итераторы
    integer nmax, n, media, i, j
    ! Массивы для полей
    real ex(ie,jb)           !поля на основной решетке
    real ey(ib,je)
    real hz(ie,je)

    real exbcf(iefbc,jebc)   !поля в передней области PML
    real eybcf(ibfbc,jebc)
    real hzxbcf(iefbc,jebc)
    real hzybcf(iefbc,jebc)

    real exbcb(iefbc,jbbc)   !поля в задней области PML
    real eybcb(ibfbc,jebc)
    real hzxbcb(iefbc,jebc)
    real hzybcb(iefbc,jebc)

    real exbcl(iebc,jb)      !поля в левой области PML
    real eybcl(iebc,je)
    real hzxbcl(iebc,je)
    real hzybcl(iebc,je)

    real exbcr(iebc,jb)      !поля в правой области PML
    real eybcr(ibbc,je)
    real hzxbcr(iebc,je)
    real hzybcr(iebc,je)

    ! Массивы для основной решетки
    real caex(ie, jb)     
    real cbex(ie, jb)
    real caey(ib, je)
    real cbey(ib, je)
    real dahz(ie, je)
    real dbhz(ie, je)

    ! Параметры областей PML
    real delbc, sigmam, bcfactor
    ! Передняя область
    real y1, y2, sigmay, ca1, cb1, caexbcf(iefbc, jebc), cbexbcf(iefbc, jebc), caexbcl(iebc, jb), &
    cbexbcl(iebc, jb), caexbcr(iebc, jb), cbexbcr(iebc, jb), sigmays, da1, db1, dahzybcf(iefbc, jebc), &
    dbhzybcf(iefbc,jebc), caeybcf(ibfbc,jebc), cbeybcf(ibfbc,jebc), dahzxbcf(iefbc,jebc), dbhzxbcf(iefbc,jebc)
    ! Задняя область
    real caexbcb(iefbc,jbbc), cbexbcb(iefbc,jbbc), dahzybcb(iefbc,jebc), dbhzybcb(iefbc,jebc), caeybcb(ibfbc,jebc), &
    cbeybcb(ibfbc,jebc), dahzxbcb(iefbc,jebc), dbhzxbcb(iefbc,jebc)
    ! Левая область
    real caeybcl(iebc, je), cbeybcl(iebc, je), x1, x2, sigmax, sigmaxs, dahzybcl(iebc, je), dbhzybcl(iebc, je), &
    dahzxbcl(iebc, je), dbhzxbcl(iebc, je)
    ! Правая область
    real caeybcr(ibbc, je), cbeybcr(ibbc, je), dahzxbcr(iebc, je), dbhzxbcr(iebc, je), dahzybcr(iebc, je), dbhzybcr(iebc, je)

    ! Динамические массивы для источника и коэффициентов обновления
    real, allocatable::source(:), ca(:), cb(:), da(:), db(:)

    ! Массивы для сохранения результатов
    real ex_res(ie, jb), ey_res(ib, je), hz_res(ie, je)
    real t_res
    ! Метка для файла сохранения и временные отсчеты
    integer results, count_start, count_end, count_rate, count_max 
    real elapsed
    character(len=32) :: arg

    ex = 0.0
    ey = 0.0
    hz = 0.0
    exbcf = 0.0
    eybcf = 0.0
    hzxbcf = 0.0
    hzybcf = 0.0
    exbcb = 0.0
    eybcb = 0.0
    hzxbcb = 0.0
    hzybcb = 0.0
    exbcl = 0.0
    eybcl = 0.0
    hzxbcl = 0.0
    hzybcl = 0.0
    exbcr = 0.0
    eybcr = 0.0
    hzxbcr = 0.0
    hzybcr = 0.0

    dx=3.0e-3
    dt=dx/(2.0*cc)  
    nmax=300 
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) nmax
    end if
    rmax=1.0e-7

    media=2
    eps=[1.0, 1.0]
    sig=[0.0, 1.0e+7]
    mur=[1.0, 1.0]
    sim=[0.0, 0.0]

    ! Излучение волны
    rtau=160.0e-12
    tau=rtau/dt
    delay=3*tau
    allocate(source(nmax))
    source=0.0
    do n=1, nmax
        source(n)=sin(omega*(n-delay)*dt)
    end do
    allocate(ca(media), cb(media), da(media), db(media))
    ! Коэффициенты обновления
    do i=1, media
        eaf=dt*sig(i)/(2.0*epsz*eps(i))
        ca(i)=(1.0-eaf)/(1.0+eaf)
        cb(i)=dt/epsz/eps(i)/dx/(1.0+eaf)
        haf=dt*sim(i)/(2.0*muz*mur(i))
        da(i)=(1.0-haf)/(1.0+haf)
        db(i)=dt/muz/mur(i)/dx/(1.0+haf)
    end do
    
    ! Геометрическое описание (основная решетка)
    ! Инициализируем всю основную решетку свободным пространством

    do i=1, ie
        do j=1, jb
            caex(i,j)=ca(1)
            cbex(i,j)=cb(1)
        end do
    end do

    do i=1, ib
        do j=1, je
            caey(i,j)=ca(1)
            cbey(i,j)=cb(1)
        end do
    end do

    do i=1, ie
        do j=1, je
            dahz(i,j)=da(1)
            dbhz(i,j)=db(1)
        end do
    end do

    ! Заполнение областей PML
    delbc=iebc*dx
    sigmam=-log(rmax)*(orderbc+1)/(2*etaz*delbc)
    bcfactor=sigmam/(dx*(orderbc+1)*(delbc**orderbc))

    ! ПЕРЕДНЯЯ область
    do i=1, iefbc
        caexbcf(i,1)=1.0
        cbexbcf(i,1)=0.0
    end do
    do j=2, jebc
        y1=(jebc-j+1.5)*dx
        y2=(jebc-j+0.5)*dx
        sigmay=bcfactor*(y1**(orderbc+1)-y2**(orderbc+1))
        ca1=exp(-sigmay*dt/epsz)
        cb1=(1.0-ca1)/(sigmay*dx)
        do i=1, iefbc
            caexbcf(i,j)=ca1
            cbexbcf(i,j)=cb1
        end do
    end do
    sigmay=bcfactor*(0.5*dx)**(orderbc+1)
    ca1=exp(-sigmay*dt/epsz)
    cb1=(1-ca1)/(sigmay*dx)
    do i=1, ie
        caex(i,1)=ca1
        cbex(i,1)=cb1
    end do
    do i=1, iebc
        caexbcl(i,1)=ca1
        cbexbcl(i,1)=cb1
    end do
    do i=1, iebc
        caexbcr(i,1)=ca1
        cbexbcr(i,1)=cb1
    end do

    do j=1, jebc
        y1=(jebc-j+1)*dx
        y2=(jebc-j)*dx
        sigmay=bcfactor*(y1**(orderbc+1)-y2**(orderbc+1))
        sigmays=sigmay*(muz/epsz)
        da1=exp(-sigmays*dt/muz)
        db1=(1-da1)/(sigmays*dx)
        do i=1, iefbc
            dahzybcf(i,j)=da1
            dbhzybcf(i,j)=db1
        end do
        do i=1, ibfbc
            caeybcf(i,j)=ca(1)
            cbeybcf(i,j)=cb(1)
        end do
        do i=1, iefbc
            dahzxbcf(i,j)=da(1)
            dbhzxbcf(i,j)=db(1)
        end do
    end do

    ! ЗАДНЯЯ область
    do i=1, iefbc
        caexbcb(i,jbbc)=1.0
        cbexbcb(i,jbbc)=0.0
    end do
    do j=2, jebc
        y1=(j-0.5)*dx
        y2=(j-1.5)*dx
        sigmay=bcfactor*(y1**(orderbc+1)-y2**(orderbc+1))
        ca1=exp(-sigmay*dt/epsz)
        cb1=(1-ca1)/(sigmay*dx)
        do i=1, iefbc
            caexbcb(i,j)=ca1
            cbexbcb(i,j)=cb1
        end do
    end do
    sigmay = bcfactor*(0.5*dx)**(orderbc+1)
    ca1=exp(-sigmay*dt/epsz)
    cb1=(1-ca1)/(sigmay*dx)
    do i=1, ie
        caex(i,jb)=ca1
        cbex(i,jb)=cb1
    end do
    do i=1, iebc
        caexbcl(i,jb)=ca1
        cbexbcl(i,jb)=cb1
    end do
    do i=1, iebc
        caexbcr(i,jb)=ca1
        cbexbcr(i,jb)=cb1
    end do

    do j=1, jebc
        y1=j*dx
        y2=(j-1)*dx
        sigmay=bcfactor*(y1**(orderbc+1)-y2**(orderbc+1))
        sigmays=sigmay*(muz/epsz)
        da1=exp(-sigmays*dt/muz)
        db1=(1-da1)/(sigmays*dx)
        do i=1, iefbc
            dahzybcb(i,j)=da1
            dbhzybcb(i,j)=db1
        end do
        do i=1, ibfbc
            caeybcb(i,j)=ca(1)
            cbeybcb(i,j)=cb(1)
        end do
        do i=1, iefbc
            dahzxbcb(i,j)=da(1)
            dbhzxbcb(i,j)=db(1)
        end do
    end do

    ! ЛЕВАЯ область 
    do j=1, je
        caeybcl(1,j)=1.0
        cbeybcl(1,j)=0.0
    end do
    do i=2, iebc
        x1=(iebc-i+1.5)*dx
        x2=(iebc-i+0.5)*dx
        sigmax=bcfactor*(x1**(orderbc+1)-x2**(orderbc+1))
        ca1=exp(-sigmax*dt/epsz)
        cb1=(1-ca1)/(sigmax*dx)
        do j=1, je
            caeybcl(i,j)=ca1
            cbeybcl(i,j)=cb1
        end do
        do j=1, jebc
            caeybcf(i,j)=ca1
            cbeybcf(i,j)=cb1
        end do
        do j=1, jebc
            caeybcb(i,j)=ca1
            cbeybcb(i,j)=cb1
        end do
    end do
    sigmax=bcfactor*(0.5*dx)**(orderbc+1)
    ca1=exp(-sigmax*dt/epsz)
    cb1=(1-ca1)/(sigmax*dx)
    do j=1, je
        caey(1,j)=ca1
        cbey(1,j)=cb1
    end do
    do j=1, jebc
        caeybcf(iebc+1,j)=ca1
        cbeybcf(iebc+1,j)=cb1
    end do
    do j=1, jebc
        caeybcb(iebc+1,j)=ca1
        cbeybcb(iebc+1,j)=cb1
    end do

    do i=1, iebc
        x1=(iebc-i+1)*dx
        x2=(iebc-i)*dx
        sigmax=bcfactor*(x1**(orderbc+1)-x2**(orderbc+1))
        sigmaxs=sigmax*(muz/epsz)
        da1=exp(-sigmaxs*dt/muz)
        db1=(1-da1)/(sigmaxs*dx)
        do j=1, je
            dahzxbcl(i,j)=da1
            dbhzxbcl(i,j)=db1
        end do
        do j=1, jebc
            dahzxbcf(i,j)=da1
            dbhzxbcf(i,j)=db1
        end do
        do j=1, jebc
            dahzxbcb(i,j)=da1
            dbhzxbcb(i,j)=db1
        end do
        do j=2, je
            caexbcl(i,j)=ca(1)
            cbexbcl(i,j)=cb(1)
        end do
        do j=1, je
            dahzybcl(i,j)=da(1)
            dbhzybcl(i,j)=db(1)
        end do
    end do

    ! ПРАВАЯ область
    do j=1, je
        caeybcr(ibbc,j)=1.0
        cbeybcr(ibbc,j)=0.0
    end do
    do i=2, iebc
        x1=(i-0.5)*dx
        x2=(i-1.5)*dx
        sigmax=bcfactor*(x1**(orderbc+1)-x2**(orderbc+1))
        ca1=exp(-sigmax*dt/epsz)
        cb1=(1-ca1)/(sigmax*dx)
        do j=1, je
            caeybcr(i,j)=ca1
            cbeybcr(i,j)=cb1
        end do
        do j=1, jebc
            caeybcf(i+iebc+ie,j)=ca1
            cbeybcf(i+iebc+ie,j)=cb1
        end do
        do j=1, jebc
            caeybcb(i+iebc+ie,j)=ca1
            cbeybcb(i+iebc+ie,j)=cb1
        end do
    end do
    sigmax=bcfactor*(0.5*dx)**(orderbc+1)
    ca1=exp(-sigmax*dt/epsz)
    cb1=(1-ca1)/(sigmax*dx)
    do j=1, je
        caey(ib,j)=ca1
        cbey(ib,j)=cb1
    end do
    do j=1, jebc
        caeybcf(iebc+ib,j)=ca1
        cbeybcf(iebc+ib,j)=cb1
    end do
    do j=1, jebc
        caeybcb(iebc+ib,j)=ca1
        cbeybcb(iebc+ib,j)=cb1
    end do

    do i=1, iebc
        x1=i*dx
        x2=(i-1)*dx
        sigmax=bcfactor*(x1**(orderbc+1)-x2**(orderbc+1))
        sigmaxs=sigmax*(muz/epsz)
        da1=exp(-sigmaxs*dt/muz)
        db1=(1-da1)/(sigmaxs*dx)
        do j=1, je
            dahzxbcr(i,j) = da1
            dbhzxbcr(i,j) = db1
        end do
        do j=1, jebc
            dahzxbcf(i+ie+iebc,j)=da1
            dbhzxbcf(i+ie+iebc,j)=db1
        end do
        do j=1, jebc
            dahzxbcb(i+ie+iebc,j)=da1
            dbhzxbcb(i+ie+iebc,j)=db1
        end do
        do j=2, je
            caexbcr(i,j)=ca(1)
            cbexbcr(i,j)=cb(1)
        end do
        do j=1, je
            dahzybcr(i,j)=da(1)
            dbhzybcr(i,j)=db(1)
        end do
    end do

    ! НАЧАЛО ЦИКЛА РАСЧЕТА РАСПРОСТРАНЕНИЯ
    call system_clock(count_start, count_rate, count_max)
    do n=1, nmax

        ! Обновление электрических полей (EX и EY) на основной сетке
        do i=1, ie
            do j=2, je
                ex(i,j)=caex(i,j)*ex(i,j)+ &
                        cbex(i,j)*(hz(i,j)-hz(i,j-1))
            end do
        end do

        do i=2, ie
            do j=1, je
                ey(i,j)=caey(i,j)*ey(i,j)+ &
                        cbey(i,j)*(hz(i-1,j)-hz(i,j))
            end do
        end do

        ! Обновление EX в областях PML
        ! ПЕРЕДНЯЯ
        do i=1, iefbc
            do j=2, jebc
                exbcf(i,j)=caexbcf(i,j)*exbcf(i,j)- &
                cbexbcf(i,j)*(hzxbcf(i,j-1)+hzybcf(i,j-1)- &
                              hzxbcf(i,j)-hzybcf(i,j))
            end do
        end do
        do i=1, ie
            ex(i,1)=caex(i,1)*ex(i,1)- &
            cbex(i,1)*(hzxbcf(ibbc+i-1,jebc)+ &
                       hzybcf(ibbc+i-1,jebc)-hz(i,1))
        end do
    
        ! ЗАДНЯЯ
        do i=1, iefbc
            do j=2, jebc-1
                exbcb(i,j)=caexbcb(i,j)*exbcb(i,j)- &
                cbexbcb(i,j)*(hzxbcb(i,j-1)+hzybcb(i,j-1)- &
                              hzxbcb(i,j)-hzybcb(i,j))
            end do
        end do
        do i=1, ie
            ex(i,jb)=caex(i,jb)*ex(i,jb)- &
            cbex(i,jb)*(hz(i,jb-1)-hzxbcb(ibbc+i-1,1)- &
                        hzybcb(ibbc+i-1,1))
        end do
        
        ! ЛЕВАЯ
        do i=1, iebc
            do j=2, je
                exbcl(i,j)=caexbcl(i,j)*exbcl(i,j)- &
                cbexbcl(i,j)*(hzxbcl(i,j-1)+hzybcl(i,j-1)- &
                              hzxbcl(i,j)-hzybcl(i,j))
            end do
        end do
        do i=1, iebc
            exbcl(i,1)=caexbcl(i,1)*exbcl(i,1)- &
            cbexbcl(i,1)*(hzxbcf(i,jebc)+hzybcf(i,jebc)- &
                          hzxbcl(i,1)-hzybcl(i,1))
        end do
        do i=1, iebc
            exbcl(i,jb)=caexbcl(i,jb)*exbcl(i,jb)- &
            cbexbcl(i,jb)*(hzxbcl(i,je)+hzybcl(i,je)- &
                           hzxbcb(i,1)-hzybcb(i,1))
        end do
        
        ! ПРАВАЯ
        do i=1, iebc
            do j=2, je
                exbcr(i,j)=caexbcr(i,j)*exbcr(i,j)- &
                cbexbcr(i,j)*(hzxbcr(i,j-1)+hzybcr(i,j-1)- &
                              hzxbcr(i,j)-hzybcr(i,j))
            end do
        end do
        do i=1, iebc
            exbcr(i,1)=caexbcr(i,1)*exbcr(i,1)- &
            cbexbcr(i,1)*(hzxbcf(iebc+ie+i,jebc)+ &
                          hzybcf(iebc+ie+i,jebc)- &
                          hzxbcr(i,1)-hzybcr(i,1))
        end do
        do i=1, iebc
            exbcr(i,jb)=caexbcr(i,jb)*exbcr(i,jb)- &
            cbexbcr(i,jb)*(hzxbcr(i,je)+hzybcr(i,je)- &
                           hzxbcb(iebc+ie+i,1)- &
                           hzybcb(iebc+ie+i,1))
        end do
        
        ! Обновление EY в областях PML

        ! ПЕРЕДНЯЯ
        do i=2, iefbc
            do j=1, jebc
                eybcf(i,j)=caeybcf(i,j)*eybcf(i,j)- &
                cbeybcf(i,j)*(hzxbcf(i,j)+hzybcf(i,j)- &
                              hzxbcf(i-1,j)-hzybcf(i-1,j))
            end do
        end do
        
        ! ЗАДНЯЯ
        do i=2, iefbc
            do j=1, jebc
                eybcb(i,j)=caeybcb(i,j)*eybcb(i,j)- &
                cbeybcb(i,j)*(hzxbcb(i,j)+hzybcb(i,j)- &
                              hzxbcb(i-1,j)-hzybcb(i-1,j))
            end do
        end do
        
        ! ЛЕВАЯ
        do i=2, iebc
            do j=1, je
                eybcl(i,j)=caeybcl(i,j)*eybcl(i,j)- &
                cbeybcl(i,j)*(hzxbcl(i,j)+hzybcl(i,j)- &
                              hzxbcl(i-1,j)-hzybcl(i-1,j))
            end do
        end do
        do j=1, je
            ey(1,j)=caey(1,j)*ey(1,j)- &
            cbey(1,j)*(hz(1,j)-hzxbcl(iebc,j)-hzybcl(iebc,j))
        end do
        
        ! ПРАВАЯ
        do i=2, iebc
            do j=1, je
                eybcr(i,j)=caeybcr(i,j)*eybcr(i,j)- &
                cbeybcr(i,j)*(hzxbcr(i,j)+hzybcr(i,j)- &
                              hzxbcr(i-1,j)-hzybcr(i-1,j))
            end do
        end do
        do j=1, je
            ey(ib,j)=caey(ib,j)*ey(ib,j)- &
            cbey(ib,j)*(hzxbcr(1,j)+hzybcr(1,j)-hz(ie,j))
        end do

        ! Обновление электрических полей (HZ) на основной сетке

        do i=1, ie
            do j=1, je
                hz(i,j)=dahz(i,j)*hz(i,j)+ &
                        dbhz(i,j)*(ex(i,j+1)-ex(i,j)+ &
                                   ey(i,j)-ey(i+1,j))
            end do
        end do
        hz(is,js)=source(n)
        ! Обновление HZX в областях PML

        ! ПЕРЕДНЯЯ
        do i=1, iefbc
            do j=1, jebc
                hzxbcf(i,j)=dahzxbcf(i,j)*hzxbcf(i,j)- &
                dbhzxbcf(i,j)*(eybcf(i+1,j)-eybcf(i,j))
            end do
        end do
        
        ! ЗАДНЯЯ
        do i=1, iefbc
            do j=1, jebc
                hzxbcb(i,j)=dahzxbcb(i,j)*hzxbcb(i,j)- &
                dbhzxbcb(i,j)*(eybcb(i+1,j)-eybcb(i,j))
            end do
        end do
        
        ! ЛЕВАЯ
        do i=1, iebc-1
            do j=1, je
                hzxbcl(i,j)=dahzxbcl(i,j)*hzxbcl(i,j)- &
                dbhzxbcl(i,j)*(eybcl(i+1,j)-eybcl(i,j))
            end do
        end do
        do j=1, je
            hzxbcl(iebc,j)=dahzxbcl(iebc,j)*hzxbcl(iebc,j)- &
            dbhzxbcl(iebc,j)*(ey(1,j)-eybcl(iebc,j))
        end do
    
        ! ПРАВАЯ
        do i=2, iebc
            do j=1, je
                hzxbcr(i,j)=dahzxbcr(i,j)*hzxbcr(i,j)- &
                dbhzxbcr(i,j)*(eybcr(i+1,j)-eybcr(i,j))
            end do
        end do
        do j=1, je
            hzxbcr(1,j)=dahzxbcr(1,j)*hzxbcr(1,j)- &
            dbhzxbcr(1,j)*(eybcr(2,j)-ey(ib,j))
        end do
        ! Обновление HZY в областях PML

        ! ПЕРЕДНЯЯ
        do i=1, iefbc
            do j=1, jebc-1
                hzybcf(i,j)=dahzybcf(i,j)*hzybcf(i,j)- &
                dbhzybcf(i,j)*(exbcf(i,j)-exbcf(i,j+1))
            end do
        end do
        do i=1, iebc
            hzybcf(i,jebc)=dahzybcf(i,jebc)*hzybcf(i,jebc)- &
            dbhzybcf(i,jebc)*(exbcf(i,jebc)-exbcl(i,1))
        end do
        do i=1, ie
            hzybcf(iebc+i,jebc)= &
            dahzybcf(iebc+i,jebc)*hzybcf(iebc+i,jebc)- &
            dbhzybcf(iebc+i,jebc)*(exbcf(iebc+i,jebc)- &
                                    ex(i,1))
        end do
        do i=1, iebc
            hzybcf(iebc+ie+i,jebc)= &
            dahzybcf(iebc+ie+i,jebc)*hzybcf(iebc+ie+i,jebc)- &
            dbhzybcf(iebc+ie+i,jebc)*(exbcf(iebc+ie+i,jebc)- &
                                       exbcr(i,1))
        end do

        ! ЗАДНЯЯ
        do i=1, iefbc
            do j=2, jebc
                hzybcb(i,j)=dahzybcb(i,j)*hzybcb(i,j)- &
                dbhzybcb(i,j)*(exbcb(i,j)-exbcb(i,j+1))
            end do
        end do
        do i=1, iebc
            hzybcb(i,1)=dahzybcb(i,1)*hzybcb(i,1)- &
            dbhzybcb(i,1)*(exbcl(i,jb)-exbcb(i,2))
        end do
        do i=1, ie
            hzybcb(iebc+i,1)= &
            dahzybcb(iebc+i,1)*hzybcb(iebc+i,1)- &
            dbhzybcb(iebc+i,1)*(ex(i,jb)-exbcb(iebc+i,2))
        end do
        do i=1, iebc
            hzybcb(iebc+ie+i,1)= &
            dahzybcb(iebc+ie+i,1)*hzybcb(iebc+ie+i,1)- &
            dbhzybcb(iebc+ie+i,1)*(exbcr(i,jb)- &
                                    exbcb(iebc+ie+i,2))
        end do
        
        ! ЛЕВАЯ
        do i=1, iebc
            do j=1, je
                hzybcl(i,j)=dahzybcl(i,j)*hzybcl(i,j)- &
                dbhzybcl(i,j)*(exbcl(i,j)-exbcl(i,j+1))
            end do
        end do
        
        ! ПРАВАЯ
        do i=1, iebc
            do j=1, je
                hzybcr(i,j)=dahzybcr(i,j)*hzybcr(i,j)- &
                dbhzybcr(i,j)*(exbcr(i,j)-exbcr(i,j+1))
            end do
        end do

        
    end do
    call system_clock(count_end, count_rate, count_max)
    elapsed = real(count_end - count_start) / real(count_rate)
    print *, 'Calculation time = ', elapsed, ' s'
    ! Сохранение полей для визуализации
    ex_res = ex
    ey_res = ey
    hz_res = hz
    t_res = dt * nmax / 1.0e-9

    open(newunit=results, file='results/fields2d.bin', &
     form='unformatted', access='stream', status='replace')
    write(results) int(ie, 4), int(jb, 4), int(ib, 4), int(je, 4), int(1, 4)
    write(results) real(dx, 4), real(dt, 4)
    write(results) real(t_res, 4)
    write(results) real(ex_res, 4)
    write(results) real(ey_res, 4)
    write(results) real(hz_res, 4)
    close(results)

end program fdtd2d