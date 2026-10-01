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
    integer, parameter::je = 300            !число ячеек сетки по оси y
    integer, parameter::ib=ie+1
    integer, parameter::jb=je+1

    integer, parameter::is=150            !расположение источника в направлении оси z
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
    real rtau, tau, delay
    real, allocatable::eaf(:), haf(:)
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
    real y1(jebc-1), y2(jebc-1), sigmay(jebc-1), ca1, cb1, ca1s(jebc-1), cb1s(jebc-1), caexbcf(iefbc, jebc), cbexbcf(iefbc, jebc), &
    caexbcl(iebc, jb), cbexbcl(iebc, jb), caexbcr(iebc, jb), cbexbcr(iebc, jb), y1e(jebc), y2e(jebc), sigmaye(jebc), sigmays(jebc), &
    da1(jebc), db1(jebc), dahzybcf(iefbc, jebc), dbhzybcf(iefbc,jebc), caeybcf(ibfbc,jebc), cbeybcf(ibfbc,jebc), dahzxbcf(iefbc,jebc), &
    dbhzxbcf(iefbc,jebc)
    ! Задняя область
    real caexbcb(iefbc,jbbc), cbexbcb(iefbc,jbbc), dahzybcb(iefbc,jebc), dbhzybcb(iefbc,jebc), caeybcb(ibfbc,jebc), &
    cbeybcb(ibfbc,jebc), dahzxbcb(iefbc,jebc), dbhzxbcb(iefbc,jebc)
    ! Левая область
    real caeybcl(iebc, je), cbeybcl(iebc, je), x1(iebc-1), x2(iebc-1), x1e(iebc), x2e(iebc), sigmax(iebc-1),sigmaxe(iebc), &
    sigmaxs(iebc), dahzybcl(iebc, je), dbhzybcl(iebc, je), dahzxbcl(iebc, je), dbhzxbcl(iebc, je)
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

    ! Массивы для замены циклов срезами
    real, allocatable::t(:)
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
    allocate(source(nmax), t(nmax))
    source=0.0  
    t = [((n - delay)*dt, n = 1, nmax)]
    source=sin(omega*t)
    allocate(ca(media), cb(media), da(media), db(media), eaf(media), haf(media))
    ! Коэффициенты обновления
    eaf=[(dt*sig(i)/(2.0*epsz*eps(i)), i = 1, media)]
    haf=[(dt*sim(i)/(2.0*muz*mur(i)), i = 1, media)]
    ca=(1.0-eaf)/(1.0+eaf)
    cb=dt/epsz/eps/dx/(1.0+eaf)
    da=(1.0-haf)/(1.0+haf)
    db=dt/muz/mur/dx/(1.0+haf)
    
    ! Геометрическое описание (основная решетка)
    ! Инициализируем всю основную решетку свободным пространством

    caex(1:ie,1:jb)=ca(1)     
    cbex(1:ie,1:jb)=cb(1)

    caey(1:ib,1:je)=ca(1)
    cbey(1:ib,1:je)=cb(1)

    dahz(1:ie,1:je)=da(1)
    dbhz(1:ie,1:je)=db(1)

    ! Заполнение областей PML
    delbc=iebc*dx
    sigmam=-log(rmax)*(orderbc+1)/(2*etaz*delbc)
    bcfactor=sigmam/(dx*(orderbc+1)*(delbc**orderbc))

    ! ПЕРЕДНЯЯ область

    caexbcf(1:iefbc,1)=1.0
    cbexbcf(1:iefbc,1)=0.0

    y1=[((jebc-j+1.5)*dx, j = 2, jebc)]
    y2=[((jebc-j+0.5)*dx, j = 2, jebc)]
    sigmay=bcfactor*(y1**(orderbc+1)-y2**(orderbc+1))
    ca1s=exp(-sigmay*dt/epsz)
    cb1s=(1.0-ca1s)/(sigmay*dx)
    caexbcf(1:iefbc,2:jebc)=spread(ca1s, 1, iefbc)
    cbexbcf(1:iefbc,2:jebc)=spread(cb1s, 1, iefbc)
    ca1=exp(-bcfactor*(0.5*dx)**(orderbc+1)*dt/epsz)
    cb1=(1-ca1)/(bcfactor*(0.5*dx)**(orderbc+1)*dx)
    caex(1:ie,1)=ca1
    cbex(1:ie,1)=cb1
    caexbcl(1:iebc,1)=ca1
    cbexbcl(1:iebc,1)=cb1
    caexbcr(1:iebc,1)=ca1
    cbexbcr(1:iebc,1)=cb1
    y1e=[((jebc-j+1)*dx, j = 1, jebc)]
    y2e=[((jebc-j)*dx, j = 1, jebc)]
    sigmaye=bcfactor*(y1e**(orderbc+1)-y2e**(orderbc+1))
    sigmays=sigmaye*(muz/epsz)
    da1=exp(-sigmays*dt/muz)
    db1=(1-da1)/(sigmays*dx)
    dahzybcf(1:iefbc,1:jebc)=spread(da1, 1, iefbc)
    dbhzybcf(1:iefbc,1:jebc)=spread(db1, 1, iefbc)
    caeybcf(1:ibfbc,1:jebc)=ca(1)
    cbeybcf(1:ibfbc,1:jebc)=cb(1)
    dahzxbcf(1:iefbc,1:jebc)=da(1)
    dbhzxbcf(1:iefbc,1:jebc)=db(1)

    ! ЗАДНЯЯ область

    caexbcb(1:iefbc,jbbc)=1.0
    cbexbcb(1:iefbc,jbbc)=0.0
    y1=[((j-0.5)*dx, j = 2, jebc)]
    y2=[((j-1.5)*dx, j = 2, jebc)]
    sigmay=bcfactor*(y1**(orderbc+1)-y2**(orderbc+1))
    ca1s=exp(-sigmay*dt/epsz)
    cb1s=(1-ca1s)/(sigmay*dx)
    caexbcb(1:iefbc,2:jebc)=spread(ca1s, 1, iefbc)
    cbexbcb(1:iefbc,2:jebc)=spread(cb1s, 1, iefbc)
    ca1=exp(-bcfactor*(0.5*dx)**(orderbc+1)*dt/epsz)
    cb1=(1-ca1)/(bcfactor*(0.5*dx)**(orderbc+1)*dx)
    caex(1:ie,jb)=ca1
    cbex(1:ie,jb)=cb1
    caexbcl(1:iebc,jb)=ca1
    cbexbcl(1:iebc,jb)=cb1
    caexbcr(1:iebc,jb)=ca1
    cbexbcr(1:iebc,jb)=cb1

    y1e=[(j*dx, j = 1, jebc)]
    y2e=[((j-1)*dx, j = 1, jebc)]
    sigmaye=bcfactor*(y1e**(orderbc+1)-y2e**(orderbc+1))
    sigmays=sigmaye*(muz/epsz)
    da1=exp(-sigmays*dt/muz)
    db1=(1-da1)/(sigmays*dx)
    dahzybcb(1:iefbc,1:jebc)=spread(da1, 1, iefbc)
    dbhzybcb(1:iefbc,1:jebc)=spread(db1, 1, iefbc)
    caeybcb(1:ibfbc,1:jebc)=ca(1)
    cbeybcb(1:ibfbc,1:jebc)=cb(1)
    dahzxbcb(1:iefbc,1:jebc)=da(1)
    dbhzxbcb(1:iefbc,1:jebc)=db(1)

    ! ЛЕВАЯ область 

    caeybcl(1,1:je)=1.0
    cbeybcl(1,1:je)=0.0
    x1=[((iebc-i+1.5)*dx, i = 2, iebc)]
    x2=[((iebc-i+0.5)*dx, i = 2, iebc)]
    sigmax=bcfactor*(x1**(orderbc+1)-x2**(orderbc+1))
    ca1s=exp(-sigmax*dt/epsz)
    cb1s=(1-ca1s)/(sigmax*dx)
    caeybcl(2:iebc,1:je)=spread(ca1s, 2, je)
    cbeybcl(2:iebc,1:je)=spread(cb1s, 2, je)
    caeybcf(2:iebc,1:jebc)=spread(ca1s, 2, jebc)
    cbeybcf(2:iebc,1:jebc)=spread(cb1s, 2, jebc)
    caeybcb(2:iebc,1:jebc)=spread(ca1s, 2, jebc)
    cbeybcb(2:iebc,1:jebc)=spread(cb1s, 2, jebc)

    ca1=exp(-bcfactor*(0.5*dx)**(orderbc+1)*dt/epsz)
    cb1=(1-ca1)/(bcfactor*(0.5*dx)**(orderbc+1)*dx)
    caey(1,1:je)=ca1
    cbey(1,1:je)=cb1
    caeybcf(iebc+1,1:jebc)=ca1
    cbeybcf(iebc+1,1:jebc)=cb1
    caeybcb(iebc+1,1:jebc)=ca1
    cbeybcb(iebc+1,1:jebc)=cb1

    x1e=[((iebc-i+1)*dx, i = 1, iebc)]
    x2e=[((iebc-i)*dx, i = 1, iebc)]
    sigmaxe=bcfactor*(x1e**(orderbc+1)-x2e**(orderbc+1))
    sigmaxs=sigmaxe*(muz/epsz)
    da1=exp(-sigmaxs*dt/muz)
    db1=(1-da1)/(sigmaxs*dx)
    dahzxbcl(1:iebc,1:je)=spread(da1, 2, je)
    dbhzxbcl(1:iebc,1:je)=spread(db1, 2, je)
    dahzxbcf(1:iebc,1:jebc)=spread(da1, 2, jebc)
    dbhzxbcf(1:iebc,1:jebc)=spread(db1, 2, jebc)
    dahzxbcb(1:iebc,1:jebc)=spread(da1, 2, jebc)
    dbhzxbcb(1:iebc,1:jebc)=spread(db1, 2, jebc)
    caexbcl(1:iebc,2:je)=ca(1)
    cbexbcl(1:iebc,2:je)=cb(1)
    dahzybcl(1:iebc,1:je)=da(1)
    dbhzybcl(1:iebc,1:je)=db(1)

    ! ПРАВАЯ область

    caeybcr(ibbc,1:je)=1.0
    cbeybcr(ibbc,1:je)=0.0
    x1=[((i-0.5)*dx, i = 2, iebc)]
    x2=[((i-1.5)*dx, i = 2, iebc)]
    sigmax=bcfactor*(x1**(orderbc+1)-x2**(orderbc+1))
    ca1s=exp(-sigmax*dt/epsz)
    cb1s=(1-ca1s)/(sigmax*dx)
    caeybcr(2:iebc,1:je)= spread(ca1s, 2, je)
    cbeybcr(2:iebc,1:je)= spread(cb1s, 2, je)
    caeybcf(ie+iebc+2 : ie+2*iebc, 1:jebc) = spread(ca1s, 2, jebc)
    cbeybcf(ie+iebc+2 : ie+2*iebc, 1:jebc) = spread(cb1s, 2, jebc)
    caeybcb(ie+iebc+2 : ie+2*iebc, 1:jebc) = spread(ca1s, 2, jebc)
    cbeybcb(ie+iebc+2 : ie+2*iebc, 1:jebc) = spread(cb1s, 2, jebc)

    ca1=exp(-bcfactor*(0.5*dx)**(orderbc+1)*dt/epsz)
    cb1=(1-ca1)/(bcfactor*(0.5*dx)**(orderbc+1)*dx)
    caey(ib,1:je)=ca1
    cbey(ib,1:je)=cb1
    caeybcf(iebc+ib,1:jebc)=ca1
    cbeybcf(iebc+ib,1:jebc)=cb1
    caeybcb(iebc+ib,1:jebc)=ca1
    cbeybcb(iebc+ib,1:jebc)=cb1

    x1e=[(i*dx, i = 1, iebc)]
    x2e=[((i-1)*dx, i = 1, iebc)]
    sigmaxe=bcfactor*(x1e**(orderbc+1)-x2e**(orderbc+1))
    sigmaxs=sigmaxe*(muz/epsz)
    da1=exp(-sigmaxs*dt/muz)
    db1=(1-da1)/(sigmaxs*dx)
    dahzxbcr(1:iebc,1:je) = spread(da1, 2, je)
    dbhzxbcr(1:iebc,1:je) = spread(db1, 2, je)
    dahzxbcf(1+ie+iebc:iefbc,1:jebc)=spread(da1, 2, jebc)
    dbhzxbcf(1+ie+iebc:iefbc,1:jebc)=spread(db1, 2, jebc)
    dahzxbcb(1+ie+iebc:iefbc,1:jebc)=spread(da1, 2, jebc)
    dbhzxbcb(1+ie+iebc:iefbc,1:jebc)=spread(db1, 2, jebc)
    caexbcr(1:iebc,2:je)=ca(1)
    cbexbcr(1:iebc,2:je)=cb(1)
    dahzybcr(1:iebc,1:je)=da(1)
    dbhzybcr(1:iebc,1:je)=db(1)

    ! НАЧАЛО ЦИКЛА РАСЧЕТА РАСПРОСТРАНЕНИЯ
    call system_clock(count_start, count_rate, count_max)
    do n=1, nmax
        ! Обновление электрических полей (EX и EY) на основной сетке
        ex(:,2:je)=caex(:,2:je)*ex(:,2:je)+ &
                cbex(:,2:je)*(hz(:,2:je)-hz(:,1:je-1))

        ey(2:ie,:)=caey(2:ie,:)*ey(2:ie,:)+ &
                cbey(2:ie,:)*(hz(1:ie-1,:)-hz(2:ie,:))

        ! Обновление EX в областях PML
        ! ПЕРЕДНЯЯ
        exbcf(:,2:jebc)=caexbcf(:,2:jebc)*exbcf(:,2:jebc)- &  
        cbexbcf(:,2:jebc)*(hzxbcf(:,1:jebc-1)+hzybcf(:,1:jebc-1)- &
                            hzxbcf(:,2:jebc)-hzybcf(:,2:jebc))
        ex(1:ie,1)=caex(1:ie,1)*ex(1:ie,1)- &
        cbex(1:ie,1)*(hzxbcf(ibbc:iebc+ie,jebc)+ &
                        hzybcf(ibbc:iebc+ie,jebc)-hz(1:ie,1))
    
        ! ЗАДНЯЯ
        exbcb(:,2:jebc-1)=caexbcb(:,2:jebc-1)*exbcb(:,2:jebc-1)- &
        cbexbcb(:,2:jebc-1)*(hzxbcb(:,1:jebc-2)+hzybcb(:,1:jebc-2)- &
                                hzxbcb(:,2:jebc-1)-hzybcb(:,2:jebc-1))
        ex(1:ie,jb)=caex(1:ie,jb)*ex(1:ie,jb)- &
        cbex(1:ie,jb)*(hz(1:ie,jb-1)-hzxbcb(ibbc:iebc+ie,1)- &
                        hzybcb(ibbc:iebc+ie,1))
        
        ! ЛЕВАЯ
        exbcl(:,2:je)=caexbcl(:,2:je)*exbcl(:,2:je)- &
        cbexbcl(:,2:je)*(hzxbcl(:,1:je-1)+hzybcl(:,1:je-1)- &
                            hzxbcl(:,2:je)-hzybcl(:,2:je))
        exbcl(:,1)=caexbcl(:,1)*exbcl(:,1)- &
        cbexbcl(:,1)*(hzxbcf(1:iebc,jebc)+hzybcf(1:iebc,jebc)- &
                        hzxbcl(:,1)-hzybcl(:,1))
        exbcl(:,jb)=caexbcl(:,jb)*exbcl(:,jb)- &
        cbexbcl(:,jb)*(hzxbcl(:,je)+hzybcl(:,je)- &
                        hzxbcb(1:iebc,1)-hzybcb(1:iebc,1))
        
        ! ПРАВАЯ
        exbcr(:,2:je)=caexbcr(:,2:je)*exbcr(:,2:je)- &
        cbexbcr(:,2:je)*(hzxbcr(:,1:je-1)+hzybcr(:,1:je-1)- &
                            hzxbcr(:,2:je)-hzybcr(:,2:je))
        exbcr(:,1)=caexbcr(:,1)*exbcr(:,1)- &
        cbexbcr(:,1)*(hzxbcf(1+iebc+ie:iefbc,jebc)+ &
                        hzybcf(1+iebc+ie:iefbc,jebc)- &
                        hzxbcr(:,1)-hzybcr(:,1))
        exbcr(:,jb)=caexbcr(:,jb)*exbcr(:,jb)- &
        cbexbcr(:,jb)*(hzxbcr(:,je)+hzybcr(:,je)- &
                        hzxbcb(1+iebc+ie:iefbc,1)- &
                        hzybcb(1+iebc+ie:iefbc,1))
        
        ! Обновление EY в областях PML

        ! ПЕРЕДНЯЯ
        eybcf(2:iefbc,:)=caeybcf(2:iefbc,:)*eybcf(2:iefbc,:)- &
        cbeybcf(2:iefbc,:)*(hzxbcf(2:iefbc,:)+hzybcf(2:iefbc,:)- &
                            hzxbcf(1:iefbc-1,:)-hzybcf(1:iefbc-1,:))
        
        ! ЗАДНЯЯ
        eybcb(2:iefbc,:)=caeybcb(2:iefbc,:)*eybcb(2:iefbc,:)- &
        cbeybcb(2:iefbc,:)*(hzxbcb(2:iefbc,:)+hzybcb(2:iefbc,:)- &
                            hzxbcb(1:iefbc-1,:)-hzybcb(1:iefbc-1,:))
        
        ! ЛЕВАЯ
        eybcl(2:iebc,:)=caeybcl(2:iebc,:)*eybcl(2:iebc,:)- &
        cbeybcl(2:iebc,:)*(hzxbcl(2:iebc,:)+hzybcl(2:iebc,:)- &
                            hzxbcl(1:iebc-1,:)-hzybcl(1:iebc-1,:))
        ey(1,:)=caey(1,:)*ey(1,:)- &
        cbey(1,:)*(hz(1,:)-hzxbcl(iebc,:)-hzybcl(iebc,:))
        
        ! ПРАВАЯ
        eybcr(2:iebc,:)=caeybcr(2:iebc,:)*eybcr(2:iebc,:)- &
        cbeybcr(2:iebc,:)*(hzxbcr(2:iebc,:)+hzybcr(2:iebc,:)- &
                            hzxbcr(1:iebc-1,:)-hzybcr(1:iebc-1,:))
        ey(ib,:)=caey(ib,:)*ey(ib,:)- &
        cbey(ib,:)*(hzxbcr(1,:)+hzybcr(1,:)- hz(ie,:))

        ! Обновление электрических полей (HZ) на основной сетке

        hz(1:ie,1:je)=dahz(1:ie,1:je)*hz(1:ie,1:je)+ & 
                    dbhz(1:ie,1:je)*(ex(1:ie,2:jb)-ex(1:ie,1:je)+ &
                                        ey(1:ie,1:je)-ey(2:ib,1:je))
        hz(is,js)=source(n)
        ! Обновление HZX в областях PML

        ! ПЕРЕДНЯЯ
        hzxbcf(1:iefbc,:)=dahzxbcf(1:iefbc,:)*hzxbcf(1:iefbc,:)- &
        dbhzxbcf(1:iefbc,:)*(eybcf(2:ibfbc,:)-eybcf(1:iefbc,:))
        
        ! ЗАДНЯЯ
        hzxbcb(1:iefbc,:)=dahzxbcb(1:iefbc,:)*hzxbcb(1:iefbc,:)- &
        dbhzxbcb(1:iefbc,:)*(eybcb(2:ibfbc,:)-eybcb(1:iefbc,:))
        
        ! ЛЕВАЯ
        hzxbcl(1:iebc-1,:)=dahzxbcl(1:iebc-1,:)*hzxbcl(1:iebc-1,:)- &
        dbhzxbcl(1:iebc-1,:)*(eybcl(2:iebc,:)-eybcl(1:iebc-1,:))
        hzxbcl(iebc,:)=dahzxbcl(iebc,:)*hzxbcl(iebc,:)- &
        dbhzxbcl(iebc,:)*(ey(1,:)-eybcl(iebc,:))
    
        ! ПРАВАЯ
        hzxbcr(2:iebc,:)=dahzxbcr(2:iebc,:)*hzxbcr(2:iebc,:)-&
        dbhzxbcr(2:iebc,:)*(eybcr(3:ibbc,:)-eybcr(2:iebc,:))
        hzxbcr(1,:)=dahzxbcr(1,:)*hzxbcr(1,:)-&
        dbhzxbcr(1,:)*(eybcr(2,:)-ey(ib,:))
        ! Обновление HZY в областях PML

        ! ПЕРЕДНЯЯ
        hzybcf(:,1:jebc-1)=dahzybcf(:,1:jebc-1)*hzybcf(:,1:jebc-1)- &
        dbhzybcf(:,1:jebc-1)*(exbcf(:,1:jebc-1)-exbcf(:,2:jebc))
        hzybcf(1:iebc,jebc)=dahzybcf(1:iebc,jebc)*hzybcf(1:iebc,jebc)- &
        dbhzybcf(1:iebc,jebc)*(exbcf(1:iebc,jebc)-exbcl(1:iebc,1))
        hzybcf(iebc+1:iebc+ie,jebc)= &
        dahzybcf(iebc+1:iebc+ie,jebc)*hzybcf(iebc+1:iebc+ie,jebc)- &
        dbhzybcf(iebc+1:iebc+ie,jebc)*(exbcf(iebc+1:iebc+ie,jebc)- &
                                        ex(1:ie,1))
        hzybcf(iebc+ie+1:iefbc,jebc)= &
        dahzybcf(iebc+ie+1:iefbc,jebc)*hzybcf(iebc+ie+1:iefbc,jebc)- &
        dbhzybcf(iebc+ie+1:iefbc,jebc)*(exbcf(iebc+ie+1:iefbc,jebc)- &
                                        exbcr(1:iebc,1))

        ! ЗАДНЯЯ
        hzybcb(1:iefbc,2:jebc)=dahzybcb(1:iefbc,2:jebc)*hzybcb(1:iefbc,2:jebc)- &
        dbhzybcb(1:iefbc,2:jebc)*(exbcb(1:iefbc,2:jebc)-exbcb(1:iefbc,3:jbbc))
        hzybcb(1:iebc,1)=dahzybcb(1:iebc,1)*hzybcb(1:iebc,1)- &
        dbhzybcb(1:iebc,1)*(exbcl(1:iebc,jb)-exbcb(1:iebc,2))
        hzybcb(iebc+1:iebc+ie,1)= &
        dahzybcb(iebc+1:iebc+ie,1)*hzybcb(iebc+1:iebc+ie,1)- &
        dbhzybcb(iebc+1:iebc+ie,1)*(ex(1:ie,jb)-exbcb(iebc+1:iebc+ie,2))
        hzybcb(iebc+ie+1:iefbc,1)= &
        dahzybcb(iebc+ie+1:iefbc,1)*hzybcb(iebc+ie+1:iefbc,1)- &
        dbhzybcb(iebc+ie+1:iefbc,1)*(exbcr(1:iebc,jb)- &
                                        exbcb(iebc+ie+1:iefbc,2))
        ! ЛЕВАЯ
        hzybcl(:,1:je)=dahzybcl(:,1:je)*hzybcl(:,1:je)- &
        dbhzybcl(:,1:je)*(exbcl(:,1:je)-exbcl(:,2:jb))
        
        ! ПРАВАЯ
        hzybcr(:,1:je)=dahzybcr(:,1:je)*hzybcr(:,1:je)- &
        dbhzybcr(:,1:je)*(exbcr(:,1:je)-exbcr(:,2:jb))
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
