    !  Srem_01.f90 
    !
    !  FUNCTIONS:
    !  Srem_01 - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Srem_01
    !
    !  PURPOSE: - Use Reference and Simulation data for the potential energy functions.
    !             Randomly select, uniformly, select a previous value.
    !
    !****************************************************************************

    module serial
    implicit none

    integer :: ios,sys_return,errnum,cmd_count
    character(120) :: cmd_arg,tmp_arg 
    integer :: run,dblevel,rep_step,temp_range,starting_temp
    integer :: InitSeed,NamdSeed
    integer :: RunLength,NamdLength,OutputFrequency,Restart,ExchangeSize
    integer :: Init_reps,InitOutputFrequency,InitSamples
    integer :: InitStartTemp,InitEndTemp,RestartFrequency,RestartStep
    double precision :: Temp_Max,Temp_Min
    integer, dimension(:), allocatable :: temperature
    double precision, dimension(:), allocatable :: temperatured
    character(len=40) :: TemplateName,FinalFile,RestartName,InitFileName
    character(len=40) :: BinaryInit,RestartNameBinary
    character(len=40) :: dummyprt
    character(len=70) :: NamdExec
    logical :: InitEnergyFromFiles,PEMethod
    character(len=40) :: PathEnergyFiles,EnergyTail,Initialize
    integer :: current_rep,current_temp,next_temp,number_atoms
    double precision :: pe_pv

    character(len=200), dimension(20) :: sedscr
    character(len=200):: exec_cmd
    character(len=316):: namd_energy_buffer

    integer, dimension(:), allocatable :: ref_pe_idx
    double precision, dimension(:,:), allocatable :: ref_pe_list
    double precision, dimension(:), allocatable :: nrxx,nryy,nrsg

    integer,parameter :: Energy_History = 1000
    integer,parameter :: Energy_Skips = 501
    character(8):: l_date
    character(10):: l_time
    double precision, parameter :: RT=0.636d0, PN=1.458397d-05
    ! PN  kcal/(mol*A^3) Chris - abeta
    ! P=2.42d-29 kcal/A^3 (1 atm)
    ! N=6.02d+23

    character*200, dimension(0:4) :: save_name,save_nameb
    integer :: save_idx = 0, exchange_type = 0, exchange_dir = 0

    double precision, dimension(:), allocatable :: relative_exch_rate
    double precision :: total_prob_of_exch
    double precision, parameter :: sqrt2 = sqrt(2.0d0)
    double precision, dimension(15) :: raw_prob_of_exchange = 0.0d0
    double precision, dimension(:), allocatable :: csum_prob_of_exchange

    type temp_array
        integer            :: index
        double precision, dimension(Energy_History)   :: temp      
    end type temp_array

    type str_temp
        integer            :: count
        double precision   :: t
        double precision   :: t2      
    end type str_temp

    type str_stat
        double precision   :: avg
        double precision   :: var
    end type str_stat

    type raw_stat
        double precision   :: sumx
        double precision   :: sumx2
    end type raw_stat

    type(str_temp), dimension(:), allocatable :: avg_temp
    type(str_stat), dimension(:), allocatable :: temp_stat
    type(temp_array), dimension(:), allocatable :: energy
    type(raw_stat), dimension(:), allocatable :: energy_stat

    contains

    ! ===================================================================================
    ! Read model parameters.
    subroutine read_parms()
    implicit none

    character*1 :: str1

    open ( unit=20, file="model_parms.txt", status='old', iostat=ios)

    read(unit=20,fmt='(i9)') dblevel

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(i7)') temp_range
    read(unit=20,fmt='(f10.5)') Temp_Min
    read(unit=20,fmt='(f10.5)') Temp_Max
    read(unit=20,fmt='(i7)') number_atoms

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(i3)') run
    read(unit=20,fmt='(i9)') rep_step
    read(unit=20,fmt='(i7)') starting_temp
    read(unit=20,fmt='(i15)') RunLength    
    read(unit=20,fmt='(i15)') NamdLength  
    read(unit=20,fmt='(i7)') OutputFrequency    
    read(unit=20,fmt='(i7)') ExchangeSize    
    read(unit=20,fmt='(l7)') PEMethod

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(a40)') TemplateName
    read(unit=20,fmt='(a40)') FinalFile 
    read(unit=20,fmt='(a40)') BinaryInit 

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(a40)') PathEnergyFiles
    read(unit=20,fmt='(a40)') EnergyTail
    read(unit=20,fmt='(a40)') Initialize

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(i7)') NamdSeed
    read(unit=20,fmt='(i7)') InitSeed

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(i7)') Restart
    read(unit=20,fmt='(i7)') RestartFrequency
    read(unit=20,fmt='(i7)') RestartStep
    read(unit=20,fmt='(a40)') RestartName
    read(unit=20,fmt='(a40)') RestartNameBinary

    read(unit=20,fmt='(a1)') str1
    read(unit=20,fmt='(a70)') NamdExec 

    close (20)

    if (dblevel .gt. 5) call log_time(1, "Parameters read in.")

    ! The random seeds must be different for each independent run.
    NamdSeed = NamdSeed * run
    InitSeed = -1 * iabs(InitSeed) * run

    ! Problem with Restart: The array in the RAN2 gets reset.  
    ! To properly fix this, the SAVE variables should be saved each frame and 
    ! reloaded upon a Restart.
    ! Quick and Dirty: Just offset the seed a little.
    if ( Restart .ne. 0 ) then
        InitSeed = InitSeed * rep_step
    endif

    end subroutine read_parms

    ! ===================================================================================
    ! Allocate arrays.
    subroutine init_arrays()
    implicit none

    integer :: ii

    allocate( temperature(temp_range) )
    allocate( temperatured(temp_range) )
    allocate( avg_temp(temp_range) )
    allocate( temp_stat(temp_range) )
    allocate( ref_pe_idx(temp_range) )
    allocate( ref_pe_list(temp_range,InitSamples+RunLength) )
    allocate( energy(temp_range) )
    allocate( energy_stat(temp_range) )
    allocate( relative_exch_rate(temp_range) )
    allocate( csum_prob_of_exchange(ExchangeSize) )

    allocate( nrxx(InitSamples) )
    allocate( nryy(InitSamples) )
    allocate( nrsg(InitSamples) )

    ! Initialize data and structures. *** Sloppy - This should be in a separate subroutine.
    temperature = 0
    temperatured = 0.0d0
    ref_pe_idx = 0
    ref_pe_list = 0.0

    do ii = 1, temp_range
        avg_temp(ii)%count = 0
        avg_temp(ii)%t = 0.0
        avg_temp(ii)%t2 = 0.0

        temp_stat(ii)%avg = 0.0
        temp_stat(ii)%var = 0.0

        energy(ii)%index = 1
        energy(ii)%temp = 0.0d0 

        energy_stat(ii)%sumx  = 0.0d0
        energy_stat(ii)%sumx2 = 0.0d0
    enddo

    nrsg = 1.0d0
    do ii = 1,InitSamples
        nrxx(ii) = dble(ii)
    enddo    

    if (dblevel .gt. 5) call log_time(2, "Arrays allocated and initialized.")

    end subroutine init_arrays

    ! ===================================================================================
    ! Allocate arrays.
    subroutine kill_arrays()

    deallocate( ref_pe_list )
    deallocate( ref_pe_idx )
    deallocate( temp_stat )
    deallocate( avg_temp )
    deallocate( temperature )    
    deallocate( temperatured )
    deallocate( nrxx )
    deallocate( nryy )
    deallocate( nrsg )
    deallocate( energy )
    deallocate( energy_stat )
    deallocate( relative_exch_rate )

    if (dblevel .gt. 5) call log_time(99, "Arrays deallocated.")

    end subroutine kill_arrays

    ! ===================================================================================
    ! Build the temperature array using exponential scaling.
    subroutine scale_temps
    implicit none

    integer :: ii
    double precision :: t1,t2

    t1 = Temp_Max * (Temp_Min/Temp_Max)**(dble(temp_range) / (dble(temp_range) - 1.0d0))
    t2 = Temp_Max / Temp_Min

    do ii = 1, temp_range
        temperatured(ii) = t1 * t2**(dble(ii)/(dble(temp_range) - 1.0d0))
        temperature(ii) = anint(t1 * t2**(dble(ii)/(dble(temp_range) - 1.0d0)))
        write( unit=98, fmt='(i3,i7)') ii,temperature(ii)
    enddo

    if (dblevel .gt. 5) call log_time(3, "Temperature array constructed.")

    end subroutine scale_temps

    ! ===================================================================================
    ! Initilalize the energy values. 
    ! *** Must initialize from files - Namd init not working.
    subroutine init_energy
    USE IFPORT

    implicit none

    integer :: ii,jj
    character*3 :: str2
    character*3 :: str3
    character*8 :: str8
    character*10 :: dummy
    character*200 :: tmpfilePE,dname
    double precision :: l_pe,l_pres,l_vol,psum

    ! Read in the initial energy distributions.
    if ( Restart .eq. 0 ) then
        if (dblevel .gt. 5) call log_time(4, "40 file read started.")

        do ii = 1, temp_range
            write (str3, '(I3.3)') temperature(ii)
            tmpfilePE = trim(PathEnergyFiles)//str3//trim(EnergyTail)

            open ( unit=30, file=tmpfilePE, status='old', iostat=ios)
            if ( ios .ne. 0 ) then
                write(unit=99,fmt='(a25,2i7)') "Open Energy File failed. ",ii,temperature(ii)
                stop
            endif
            do jj = 1, Energy_Skips
                read( unit=30, fmt='(a1)', iostat=ios ) dummy
            enddo

            do jj = 1, Energy_History
                read( unit=30, fmt='(190x,f15.4,35x,f15.4,15x,f15.4)', iostat=ios ) l_pe,l_pres,l_vol
                energy(ii)%temp(jj) = l_pe + l_vol*PN
            enddo
            close (30) 
        enddo

        if (dblevel .gt. 5) call log_time(5, "40 file read ended for restart.")

    else
        ! Read the binary energy data for the run.
        write (str2, '(I3.3)') run
        write (str8, '(I8.8)') RestartStep
        tmpfilePE = "./output/Energy_"//trim(FinalFile)//"_r"//str2//str8//".bin"

        open ( unit=31, file=tmpfilePE, status='old',form='binary', iostat=ios)
        if ( ios .ne. 0 ) then
            write(unit=99,fmt='(a36)') "Open Energy File failed for restart."
            dname = adjustl(tmpfilePE)
            write(unit=99,fmt='(a200)') dname
            stop
        endif
        read(unit=31) energy
        close (31)      
    endif

    if (dblevel .gt. 5) call log_time(6, "Energy data read in.")

    ! Data for the running mean and variances of the energy.
    do ii = 1, temp_range
        energy_stat(ii)%sumx  = sum( energy(ii)%temp(:) )
        energy_stat(ii)%sumx2 = sum( energy(ii)%temp(:)**2 )
    enddo

    call write_energy_stats(0)

    ! Get the probability of exchanging to multiple temperatures.
    tmpfilePE = trim(Initialize)//"exchange_prob.txt"

    open ( unit=32, file=tmpfilePE, status='old', iostat=ios)
    if ( ios .ne. 0 ) then
        write(unit=99,fmt='(a39)') "Openning exchange_prob.txt file failed."
        dname = adjustl(tmpfilePE)
        write(unit=99,fmt='(a200)') dname
        stop
    endif

    do ii = 1,15
        read(unit=32,fmt='(f12.8)', iostat=ios) raw_prob_of_exchange(ii)
        if ( ios .ne. 0 ) then
            write(unit=99,fmt='(a52,i5)') "Problem reading the exchange_prob.txt file at line: ", ii
            stop
        endif
    enddo
    close (32) 

    ! Normalize the probabilities over the maximum size of the exchanges.
    psum = 0.0d0
    do ii = 1, ExchangeSize
        psum = psum + raw_prob_of_exchange(ii)
    enddo

    ! The cummulative sum will be compared to a random uniform to find the jump size.
    csum_prob_of_exchange = 0.0d0
    csum_prob_of_exchange(1) = raw_prob_of_exchange(1) / psum
    do ii = 2, ExchangeSize
        csum_prob_of_exchange(ii) = csum_prob_of_exchange(ii-1) + raw_prob_of_exchange(ii) / psum
    enddo

    do ii = 1, ExchangeSize
        write( unit=98, fmt='(i5,2f12.8)') ii,raw_prob_of_exchange(ii),csum_prob_of_exchange(ii)
    enddo

    if (dblevel .gt. 5) call log_time(7, "Probabilities for multiple exchanges.")

    end subroutine init_energy

    ! ===================================================================================
    ! Run the replicas.
    subroutine run_reps()

    USE IFPORT

    implicit none

    integer :: rep,exch=1
    integer :: rep_seed
    integer :: vidx,record_count,run_factor
    double precision :: scale
    character*8 :: str1
    character*3 :: str2
    character*3 :: str7
    character*10 :: str3,str4,str5,str6
    character*200 :: syscall,tmpfileIn,tmpfileOut
    character*200 :: rep_name,dcd_name,old_rep_name
    character*200 :: rep_nameb,old_rep_nameb,leftalign
    character*80 :: lineFull
    character*4 :: lineID
    character*30 :: lineA
    character*26 :: lineB
    double precision :: velx,vely,velz

    current_temp = starting_temp

    !if ( rep_step .eq. 1 ) then
    !    syscall = "rm energy*.dat"
    !    sys_return = system( syscall )
    !endif

    do rep = rep_step, RunLength
        ! Each rep gets a different random number seed.
        rep_seed = NamdSeed + rep        
        current_rep = rep

        ! Convert integers to character strings.
        write (str1, '(I8.8)') rep
        write (str2, '(I3.3)') current_temp
        write (str3, '(I10)') temperature(current_temp)
        write (str4, '(I10)') NamdLength
        write (str5, '(I10)') OutputFrequency
        write (str6, '(I10)') rep_seed
        write (str7, '(I3.3)') run

        rep_name  = trim(FinalFile)//"_r"//trim(str7)//trim(str1)//"_t"//trim(str2)
        rep_nameb = trim(FinalFile)//"i_r"//trim(str7)//trim(str1)//"_t"//trim(str2)
        dcd_name  = trim(FinalFile)//"_r"//trim(str7)//trim(str1)//"_t"//trim(str2)

        open ( unit=30, file="sedscr.txt", status='replace', iostat=ios)

        write( unit=30,fmt='(a)') trim("s/FinalFile/"//trim(rep_name)//"/g")
        write( unit=30,fmt='(a)') trim("s/RepTemp/"//trim(adjustl(str3))//"/g")
        write( unit=30,fmt='(a)') trim("s/NamdLength/"//trim(adjustl(str4))//"/g")
        write( unit=30,fmt='(a)') trim("s/RestartFile/"//trim(rep_nameb)//"/g")
        write( unit=30,fmt='(a)') trim("s/OutputFrequency/"//trim(adjustl(str5))//"/g")
        write( unit=30,fmt='(a)') trim("s/RandSeed/"//trim(adjustl(str6))//"/g")
        write( unit=30,fmt='(a)') trim("s/DCDFile/"//trim(dcd_name)//"/g")

        close(30)

        ! Edit the template file.
        syscall = trim( "sed -f sedscr.txt "//trim(TemplateName)//" > "//trim(rep_name)//".namd" )
        sys_return = system( syscall )
        if ( sys_return /= 0 ) then
            errnum = ierrno( )
            write( unit=99, fmt='(a44,i7)') "Error in initialization sed call: errnum = ", errnum
            write( unit=99, fmt='(a15,i7,a20,i7)') "Init_Reps = ", rep ,"    Temp_Range = ", current_temp
            stop
        end if

        write( dummyprt, fmt='(a17,i8,a7,2i4)') "Template Rep: ",rep, " Temp: ", current_temp, temperature(current_temp)
        if (dblevel .gt. 20) call log_time(20, dummyprt)

        ! Move input files, rep_name contains the new temperature index.
        if ( rep .eq. 1 ) then                ! Start of run.
            syscall = "cp -f "//trim(Initialize)//trim(FinalFile)//trim(adjustl(str3))//".coor  ./output/InText.coor"
            sys_return = system( syscall )
            syscall = "cp -f "//trim(Initialize)//trim(BinaryInit)//trim(adjustl(str3))//".coor ./output/InBinary.coor"
            sys_return = system( syscall )
            syscall = "cp -f "//trim(Initialize)//trim(BinaryInit)//trim(adjustl(str3))//".vel  ./output/InBinary.vel"
            sys_return = system( syscall )
            syscall = "cp -f "//trim(Initialize)//trim(BinaryInit)//trim(adjustl(str3))//".xsc  ./output/InBinary.xsc"
            sys_return = system( syscall )
        else if ( Restart .eq. 0 ) then       ! Normal continuation.
            syscall = "cp -f ./output/"//trim(old_rep_name)//".coor   ./output/InText.coor"
            sys_return = system( syscall )
            syscall = "cp -f ./output/"//trim(old_rep_nameb)//".coor  ./output/InBinary.coor"
            sys_return = system( syscall )
            syscall = "cp -f ./output/"//trim(old_rep_nameb)//".vel   ./output/InBinary.vel"
            sys_return = system( syscall )
            syscall = "cp -f ./output/"//trim(old_rep_nameb)//".xsc   ./output/InBinary.xsc"
            sys_return = system( syscall )
        else ! Restart.
            Restart = 0
            syscall = "cp -f ./output/"//trim(RestartName)//".coor       ./output/InText.coor"
            sys_return = system( syscall )
            syscall = "cp -f ./output/"//trim(RestartNameBinary)//".coor ./output/InBinary.coor"
            sys_return = system( syscall )
            syscall = "cp -f ./output/"//trim(RestartNameBinary)//".vel  ./output/InBinary.vel"
            sys_return = system( syscall )
            syscall = "cp -f ./output/"//trim(RestartNameBinary)//".xsc  ./output/InBinary.xsc"
            sys_return = system( syscall )
        endif

        if (dblevel .gt. 20) call log_time(21, "Input files moved.")

        ! Run Namd.
        tmpfileIn = trim(rep_name)//".namd"
        tmpfileOut = trim(rep_name)//".out"
        syscall = trim( trim(NamdExec)//" "//trim(tmpfileIn)//" > "//trim(tmpfileOut) )
        sys_return = system( syscall )

        ! If namd failed, try running it again.
        if ( sys_return /= 0 ) then
            errnum = ierrno( )
            syscall = "mv -f "//trim(tmpfileOut)//"  "//trim(tmpfileOut)//".BAK"
            sys_return = system( syscall )
            write( unit=99, fmt='(a33,i7)') "Error in running Namd: errnum = ", sys_return
            write( unit=99, fmt='(a15,i7,a20,i7)') "Init_Reps = ", rep ,"    Temp_Range = ", current_temp
            leftalign = adjustl(trim(syscall))
            write( unit=99, fmt='(a200)') leftalign
            syscall = trim( trim(NamdExec)//" "//trim(tmpfileIn)//" > "//trim(tmpfileOut) )
            sys_return = system( syscall )
        end if

        ! If namd fails a second, try running it again.
        if ( sys_return /= 0 ) then
            errnum = ierrno( )
            syscall = "mv -f "//trim(tmpfileOut)//"  "//trim(tmpfileOut)//".BAK2"
            sys_return = system( syscall )
            write( unit=99, fmt='(a33,i7)') "Error in running Namd: errnum = ", sys_return
            write( unit=99, fmt='(a15,i7,a20,i7)') "Init_Reps = ", rep ,"    Temp_Range = ", current_temp
            leftalign = adjustl(trim(syscall))
            write( unit=99, fmt='(a200)') leftalign
            syscall = trim( trim(NamdExec)//" "//trim(tmpfileIn)//" > "//trim(tmpfileOut) )
            sys_return = system( syscall )
        end if

        if ( sys_return /= 0 ) then
            errnum = ierrno( )
            write( unit=99, fmt='(a44,i7)') "Error in initialization Namd call: errnum = ", sys_return
            write( unit=99, fmt='(a15,i7,a20,i7)') "Init_Reps = ", rep ,"    Temp_Range = ", current_temp
            stop
        end if          

        write( dummyprt, fmt='(a17,i8,a7,2i4)') "Namd exits. Rep: ",rep, " Temp: ", current_temp, temperature(current_temp)
        if (dblevel .gt. 20) call log_time(22, dummyprt)

        ! Bookkeeping.
        old_rep_name  = trim(rep_name)
        old_rep_nameb = trim(rep_nameb)

        ! Process the Namd output.
        call get_energy(tmpfileOut)

        ! Check for exchange.
        exch = exchange()

        ! Update the energy array.
        call update_energy_array()

        ! Log the type of exchange.
        write( dummyprt, fmt='(a6,i4,i8,4i5)') "Exch: ",run,rep,current_temp,exch,exchange_type,exchange_dir
        if (dblevel .gt. 15) call log_time(23, dummyprt)

        if ( exch .ne. 0 ) then
            ! Scale the velocities to the new temperature.
            next_temp = current_temp + exch

            scale = sqrt( dble(temperature(next_temp)) / dble(temperature(current_temp)) )

            ! Only the one velocity file needs modification.
            tmpfileIn  = "./output/"//trim(old_rep_nameb)//".vel"
            tmpfileOut = "./output/tempexch.vel"

            open ( unit=50, file=trim(tmpfileIn), status='old',form='binary', iostat=ios)
            open ( unit=51, file=trim(tmpfileOut), status='replace',form='binary', iostat=ios)

            read ( unit=50 ) record_count
            if ( record_count .ne. number_atoms ) then
                write ( unit=99,fmt='(a47)') "Incorrect number of atoms in the velocity file."
                write ( unit=99,fmt='(a10,i8,a11,i8)') "Expected: ", number_atoms,"  In file: ", record_count
                write( unit=99, fmt='(a15,i7,a20,i7)') "Reps = ", rep ,"    Temp = ", current_temp
                stop
            endif

            write ( unit= 51 ) record_count

            do vidx = 1,record_count
                read ( unit=50 ) velx,vely,velz
                write( unit=51 ) velx*scale,vely*scale,velz*scale
            enddo

            close (50)
            close (51)

            syscall = "rm -f "//trim(tmpfileIn)
            sys_return = system( syscall )
            syscall = "mv -f "//trim(tmpfileOut)//" "//trim(tmpfileIn)
            sys_return = system( syscall )

        endif

        current_temp = current_temp + exch

        call cleanup(rep,old_rep_name,old_rep_nameb)

        if ( mod(rep,100) .eq. 0 ) call write_energy_stats(rep)

    enddo

    if (dblevel .gt. 25) call log_time(98, "Replication loop complete.")

    end subroutine run_reps

    ! ===================================================================================
    ! Return an integer for the direction and size of the exchange.
    ! Exchange_type reports where in the algorithm the return occurs.
    Function exchange()
    implicit none

    integer :: exchange
    integer :: direction,rand_idx
    double precision, parameter :: kb = 0.0019872041d0
    double precision :: potential, pressure, volume
    double precision :: randenergy, delta, expdelta
    double precision :: rdir,rnorm,rmet

    exchange = 0
    exchange_type = 0

    read(namd_energy_buffer,fmt='(190x,f15.4,35x,f15.4,15x,f15.4)') potential, pressure, volume

    pe_pv = potential + volume*PN    

    ! Is the move up or down?
    rdir = ran2(InitSeed)
    if ( rdir .GE. 0.5d0 ) then
        direction = 1
    else
        direction = -1
    endif

    exchange_type = 1

    if ( ExchangeSize .gt. 1 ) then
        call Multiple_Temps(direction)
    endif

    exchange_dir = direction

    ! No moves beyond boundaries.
    if ( current_temp + direction .lt. 1 .or. &
    current_temp + direction .gt. temp_range ) return

    exchange_type = 2

    ! Pick a random energy at the target temperature.
    !rnorm = gasdev(InitSeed)
    !randenergy = temp_stat(current_temp + direction).avg + rnorm * sqrt(temp_stat(current_temp + direction).var)

    rand_idx = int( ran2(InitSeed) * 1000 ) + 1
    randenergy = energy(current_temp + direction)%temp(rand_idx)    

    delta = (pe_pv - randenergy) * (1.0d0/(kb*temperature(current_temp)) - 1.0d0/(kb*temperature(current_temp + direction)))

    ! Always exchange if the energy at the higher temperature is lower than current.
    if ( delta .ge. 0.0d0 ) then
        exchange = direction
        exchange_type = 3
        return
    endif

    exchange_type = 4

    ! Continue Metropolis criteria.
    rmet = ran2(InitSeed)
    expdelta = exp(delta)
    if ( rmet .lt. expdelta ) then
        exchange = direction
        exchange_type = 5
        return
    endif

    exchange_type = 6

    if (dblevel .gt. 25) call log_time(24, "Exchange complete.")

    return

    end function exchange

    ! ===================================================================================
    ! Determine if the exchange should cross cross multiple temperatures.
    ! 
    subroutine Multiple_Temps(direction)

    integer, intent(inout) :: direction
    integer :: ii, direction_sign
    double precision :: current_mean,current_var,test_mean,test_var
    double precision :: random_selection
    double precision :: cum_sum

    ! Keep track if the move is up or down.
    direction_sign = direction

    ! Immediately return if over a boundary.
    if ( direction .gt. 0 ) then
        ! Can't move over the upper boundary.
        if ( current_temp .eq. temp_range ) return
    else
        ! Can't move over the lower boundary.
        if ( current_temp .eq. 1 ) return
    endif

    ! Randomly select a temperature for an attempted exchange.
    random_selection = ran2(InitSeed)

    do ii = 1, ExchangeSize        
        if ( random_selection .le. csum_prob_of_exchange(ii) ) then
            direction = sign (ii, direction_sign )
            return
        endif
    enddo

    ! One of the temperatures should have been selected.  Hard crash.
    write( unit=99, fmt='(a39)' ) "Error in the Multiple_Temps subroutine."
    write( unit=99, fmt='(a18,f15.12)' ) "Random_Selection: ", random_selection
    write( unit=99, fmt='(a18)' ) "relative_exch_rate"
    write( unit=99, fmt='(100f10.6)' ) relative_exch_rate
    stop

    end subroutine Multiple_Temps

    ! ===================================================================================
    ! From Numerical Recipes
    FUNCTION ran1(idum)
    INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
    REAL ran1,AM,EPS,RNMX
    PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/, iy /0/
    if (idum.le.0.or.iy.eq.0) then
        idum=max0(-idum,1)
        do j=NTAB+8,1,-1
            k=idum/IQ
            idum=IA*(idum-k*IQ)-IR*k
            if (idum.lt.0) idum=idum+IM
            if (j.le.NTAB) iv(j)=idum
        enddo
        iy=iv(1)
    endif
    k=idum/IQ
    idum=IA*(idum-k*IQ)-IR*k
    if (idum.lt.0) idum=idum+IM
    j=1+iy/NDIV
    iy=iv(j)
    iv(j)=idum
    ran1=min(AM*iy,RNMX)
    return
    END FUNCTION ran1

    ! ===================================================================================
    ! From Numerical Recipes
    FUNCTION ran2(idum)
    INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL ran2,AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
    NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
    endif
    k=idum/IQ1
    idum=IA1*(idum-k*IQ1)-k*IR1
    if (idum.lt.0) idum=idum+IM1
    k=idum2/IQ2
    idum2=IA2*(idum2-k*IQ2)-k*IR2
    if (idum2.lt.0) idum2=idum2+IM2
    j=1+iy/NDIV
    iy=iv(j)-idum2
    iv(j)=idum
    if(iy.lt.1)iy=iy+IMM1
    ran2=min(AM*iy,RNMX)
    return
    END FUNCTION ran2

    ! ===================================================================================
    ! From Numerical Recipes
    FUNCTION gasdev(idum)
    INTEGER idum
    REAL gasdev
    ! CU    USES ran1
    INTEGER iset
    REAL fac,gset,rsq,v1,v2!,ran1
    SAVE iset,gset
    DATA iset/0/
    if (idum.lt.0) iset=0
    if (iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
    else
        gasdev=gset
        iset=0
    endif
    return
    END FUNCTION gasdev

    ! ===================================================================================
    ! From Numerical Recipes
    SUBROUTINE fit(x,y,ndata,sig,mwt,a,b,siga,sigb,chi2,q)
    INTEGER :: mwt,ndata
    double precision :: a,b,chi2,q,siga,sigb,sig(ndata),x(ndata),y(ndata)
    INTEGER :: i
    double precision :: sigdat,ss,st2,sx,sxoss,sy,t,wt
    sx=0.
    sy=0.
    st2=0.
    b=0.
    if(mwt.ne.0) then
        ss=0.
        do i=1,ndata
            wt=1./(sig(i)**2)
            ss=ss+wt
            sx=sx+x(i)*wt
            sy=sy+y(i)*wt
        enddo
    else
        do i=1,ndata
            sx=sx+x(i)
            sy=sy+y(i)
        enddo
        ss=float(ndata)
    endif
    sxoss=sx/ss
    if(mwt.ne.0) then
        do i=1,ndata
            t=(x(i)-sxoss)/sig(i)
            st2=st2+t*t
            b=b+t*y(i)/sig(i)
        enddo
    else
        do i=1,ndata
            t=x(i)-sxoss
            st2=st2+t*t
            b=b+t*y(i)
        enddo
    endif
    b=b/st2
    a=(sy-sx*b)/ss
    siga=sqrt((1.+sx*sx/(ss*st2))/ss)
    sigb=sqrt(1./st2)
    chi2=0.
    q=1.
    if(mwt.eq.0) then
        do i=1,ndata
            chi2=chi2+(y(i)-a-b*x(i))**2
        enddo
        sigdat=sqrt(chi2/(ndata-2))
        siga=siga*sigdat
        sigb=sigb*sigdat
    else
        do i=1,ndata
            chi2=chi2+((y(i)-a-b*x(i))/sig(i))**2
        enddo
        if(ndata.gt.2) then
            q = gammq(0.5d0*(ndata-2),0.5*chi2)
        endif
    endif
    return
    end subroutine fit

    ! ===================================================================================
    ! From Numerical Recipes
    FUNCTION gammq(a,x)
    double precision :: a,gammq,x
    double precision :: gammcf,gamser,gln
    if(x.lt.0..or.a.le.0.) then
        write( unit=99, fmt='(a28,f10.4,a7,f10.4)') "Bad arguments in gammq: x = ",x,"   a = ", a
    endif
    if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
    else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
    endif
    return
    END FUNCTION gammq

    ! ===================================================================================
    ! From Numerical Recipes
    SUBROUTINE gcf(gammcf,a,x,gln)
    double precision ::  a,gammcf,gln,x
    INTEGER, PARAMETER :: ITMAX=100
    double precision, PARAMETER :: EPS=3.e-7,FPMIN=1.e-30
    INTEGER :: i
    double precision :: an,b,c,d,del,h
    gln=gammln(a)
    b=x+1.-a
    c=1./FPMIN
    d=1./b
    h=d
    do i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS) then 
            gammcf=exp(-x+a*log(x)-gln)*h
            return            
        endif
    enddo
    write( unit=99, fmt='(a35)') "a too large, ITMAX too small in gcf"
    stop
    end subroutine gcf

    ! ===================================================================================
    ! From Numerical Recipes
    SUBROUTINE gser(gamser,a,x,gln)
    INTEGER, PARAMETER :: ITMAX=100
    double precision :: a,gamser,gln,x
    double precision, PARAMETER :: EPS=3.e-7
    INTEGER :: n
    double precision :: ap,del,sum
    gln=gammln(a)
    if(x.le.0.)then
        if(x.lt.0.) then
            write( unit=99, fmt='(a13)') "x < 0 in gser"
        endif
        gamser=0.
        return
    endif
    ap=a
    sum=1./a
    del=sum
    do n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS) then 
            gamser=sum*exp(-x+a*log(x)-gln)
            return
        endif
    enddo

    write( unit=99, fmt='(a35)') "a too large, ITMAX too small in gser"
    stop

    end subroutine gser

    ! ===================================================================================
    ! From Numerical Recipes
    FUNCTION gammln(xx)
    double precision :: gammln,xx
    INTEGER :: j
    double precision :: ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0, &
    +24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
    -.5395239384953d-5,2.5066282746310005d0/
    x=xx
    y=x
    tmp=x+5.5d0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=1.000000000190015d0
    do j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    return
    END FUNCTION gammln

    ! ===================================================================================
    ! From Numerical Recipes
    FUNCTION erfcc(x)
    double precision :: erfcc,x
    double precision :: t,z
    z = abs(x)
    t = 1.0d0/(1.0d0 + 0.5d0 * z)
    erfcc = t * exp(-z * z-1.26551223+t*(1.00002368+t*(.37409196+t* &
    (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t* &
    (1.48851587+t*(-.82215223+t*.17087277)))))))))
    if (x .lt. 0.0) erfcc = 2.0d0 - erfcc
    return
    END FUNCTION erfcc

    ! ===================================================================================
    ! Print the time and date to the log with a message.
    subroutine log_time(nn,aa)

    integer, intent(in) :: nn
    character(len=*), intent(in) :: aa
    character(len=2) :: str2

    write (str2, '(I2.2)') nn
    call date_and_time(l_date,l_time)    
    write( unit=98, fmt='(a2,2x,a40,2x,a8,2x,a10)') str2,aa,l_date,l_time

    end subroutine log_time

    ! ===================================================================================
    ! Read the last ENERGY line from the NAMD output file.
    ! Write this to the log and global buffer.
    subroutine get_energy(tmpfileOut)

    character(len=6) :: buffer
    character*200, intent(in) :: tmpfileOut

    open ( unit=81, file=trim(tmpfileOut), status='old', iostat=ios)

    do
        read(unit=81,fmt='(a6)',end=100) buffer
        if ( buffer .ne. "ENERGY" ) cycle
        backspace(81)
        read(unit=81,fmt='(a316)') namd_energy_buffer
    enddo
100 continue    

    write(unit=97,fmt='(i3,1x,i8,1x,i3,1x,a316)') run,current_rep,current_temp,namd_energy_buffer
    flush(97)

    close(81)

    end subroutine get_energy

    ! ===================================================================================
    ! Update the energy data.  It is in a circular queue, so loop over and restart at the end.
    subroutine update_energy_array()

    integer :: next_index
    character(3) :: str3
    character(8) :: str8
    character(120) :: energyfile
    double precision :: old_energy

    ! Energy_History is the number of elements, i.e. 1000
    next_index = mod( energy(current_temp)%index, Energy_History ) + 1

    ! The pe_pv is in a global variable and was updated in exchange()

    ! First remove the last element's data from the running history.
    old_energy = energy(current_temp)%temp(energy(current_temp)%index)
    energy_stat(current_temp)%sumx  = energy_stat(current_temp)%sumx  - old_energy
    energy_stat(current_temp)%sumx2 = energy_stat(current_temp)%sumx2 - old_energy**2

    ! Then add the latest energy back in.
    energy_stat(current_temp)%sumx  = energy_stat(current_temp)%sumx  + pe_pv
    energy_stat(current_temp)%sumx2 = energy_stat(current_temp)%sumx2 + pe_pv**2

    ! The latest energy also goes into the queue.
    energy(current_temp)%temp(energy(current_temp)%index) = pe_pv

    ! Point to the next element to be updated.
    energy(current_temp)%index = next_index

    write (str3, '(I3.3)') run
    write (str8, '(I8.8)') current_rep
    energyfile = "./output/Energy_"//trim(FinalFile)//"_r"//str3//str8//".bin"

    if (dblevel .gt. 30) call log_time(30, "Binary write started.")

    open ( unit=31, file=triM(energyfile), status='replace',form='binary', iostat=ios)
    write(unit=31) energy
    close (31) 

    if (dblevel .gt. 30) call log_time(31, "Binary write ended.")

    end subroutine update_energy_array

    ! ===================================================================================
    ! Delete old coor,vel,xsc, and bin files.
    ! All dcd files are retained, and all others with multiples of the RestartFrequency.
    ! The deletions start 5 back from the current rep, so the recent runs have a 
    ! delayed deletion.  The Energy~bin files need the temperature removed and Energy_ added.
    subroutine cleanup(rep,old_rep_name,old_rep_nameb)
    USE IFPORT

    integer, intent(in) :: rep    
    character*200, intent(in) :: old_rep_name,old_rep_nameb
    character*200 :: syscall,r_name,l_name
    integer :: to_delete,sys_return

    ! This is the oldest element in the circular queue.
    to_delete = mod(save_idx+4, 5)

    ! The system must be 5 reps in, and not a multiple of the restart frequency to be deleted.
    if ( rep - rep_step .ge. 5 .and. mod(rep-5,RestartFrequency) .ne. 0) then
        syscall = "rm -f ./output/"//trim(save_name(to_delete))//".coor"
        sys_return = system( syscall )
        syscall = "rm -f ./output/"//trim(save_name(to_delete))//".vel"
        sys_return = system( syscall )
        syscall = "rm -f ./output/"//trim(save_name(to_delete))//".xsc"
        sys_return = system( syscall )
        syscall = "rm -f ./output/"//trim(save_nameb(to_delete))//".coor"
        sys_return = system( syscall )
        syscall = "rm -f ./output/"//trim(save_nameb(to_delete))//".vel"
        sys_return = system( syscall )
        syscall = "rm -f ./output/"//trim(save_nameb(to_delete))//".xsc"
        sys_return = system( syscall )        
        r_name = adjustr(save_name(to_delete))
        l_name = adjustl( r_name( 1 : 200-5 ) )
        syscall = "rm -f ./output/Energy_"//trim(l_name)//".bin"
        sys_return = system( syscall )
    endif

    ! Add the latest file names to the queue.
    save_name(to_delete) = old_rep_name
    save_nameb(to_delete) = old_rep_nameb

    ! Increment the queue index.
    save_idx = mod(save_idx+1, 5)

    if (dblevel .gt. 20) call log_time(25, "Output files deleted.")

    end subroutine cleanup

    ! ===================================================================================
    ! Write the energy statistics.
    subroutine write_energy_stats(rep)

    integer, intent(in) :: rep
    integer :: ii

    write(unit=96,fmt='(i9,100f14.5)') rep, ( energy_stat(ii)%sumx/dble(Energy_History), &
    sqrt(energy_stat(ii)%sumx2/dble(Energy_History) - (energy_stat(ii)%sumx/dble(Energy_History))**2), &
    ii = 1, temp_range )

    end subroutine write_energy_stats


    ! ===================================================================================

    end module serial

    ! ===================================================================================
    ! ===================================================================================

    program Srem_01

    use serial
    USE IFPORT

    implicit none

    ! Variables
    integer :: ii,jj
    character(80) :: fileout

    !! Command line interface.
    cmd_count = command_argument_count()
    call get_command_argument(0,cmd_arg)

    ! Standard log files.

    open ( unit=99, file="error_log.txt", status='replace', iostat=ios)
    !open ( unit=97, file="namd_energy.txt", status='replace', iostat=ios)
    open ( unit=97, file="namd_energy.txt", status='new', iostat=ios)
    if ( ios .ne. 0 ) then
        write(unit=99,fmt='(a79)') "Existing namd_energy.txt file found.  This must be manually deleted or renamed."
        stop
    endif
    open ( unit=98, file="run_log.txt", status='replace', iostat=ios)
    open ( unit=96, file="stats_energy.txt", status='replace', iostat=ios)

    ! Write the command line to the log.
    tmp_arg = adjustl(cmd_arg)
    write( unit=98,fmt='(a120)') tmp_arg

    ! Read in parameters and set up arrays.
    call read_parms()
    call init_arrays()

    ! Still initializing.
    call scale_temps()
    call init_energy()   

    ! Run namd. 
    call run_reps()

    ! Clean up.
    call kill_arrays()
    close (96)
    close (97)
    close (98)
    close (99)

    end program Srem_01


