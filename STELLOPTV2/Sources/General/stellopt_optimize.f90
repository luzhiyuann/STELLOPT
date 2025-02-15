!-----------------------------------------------------------------------
!     Subroutine:    stellopt_optimize
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          05/26/2012
!     Description:   This subroutine call the optimization routine.
!-----------------------------------------------------------------------
      SUBROUTINE stellopt_optimize
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stellopt_runtime
      USE stellopt_input_mod
      USE stellopt_vars
      USE safe_open_mod, ONLY: safe_open
      USE fdjac_mod, ONLY: FLAG_CLEANUP, FLAG_CLEANUP_LEV, FLAG_SINGLETASK
      USE mpi_params
      USE mpi_inc
      
!-----------------------------------------------------------------------
!     Local Variables
!        ier         Error flag
!        iunit       File unit number
!----------------------------------------------------------------------
      IMPLICIT NONE
      !LOGICAL ::  lrestart
      LOGICAL ::  lfile_exists, lskip_min, ldeleteopt, lsocleanup
      INTEGER ::  ier, iunit,nvar_in, nprint, info, ldfjac,nfev,&
                  iunit_restart, nfev_save, npop, ndiv, i
      INTEGER, ALLOCATABLE :: ipvt(:)
      REAL(rprec)              ::  target_fitness, c1, c2
      REAL(rprec), ALLOCATABLE ::  qtf(:), wa1(:), wa2(:), wa3(:), &
                                   wa4(:), fvec(:)
      REAL(rprec), ALLOCATABLE ::  fjac(:,:)
      LOGICAL :: used_mango_algorithm
      REAL(rprec), EXTERNAL :: enorm
      CHARACTER(200) :: cmdtxt = ""
      EXTERNAL stellopt_fcn
      
!----------------------------------------------------------------------
!     BEGIN SUBROUTINE
!----------------------------------------------------------------------
!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)                       !PPPL -SAL
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_optimize1',ierr_mpi)
!DEC$ ENDIF
      IF (npopulation <= 0) npopulation = numprocs
      CALL tolower(opt_type)

      ! Print to screen
      IF (lverb) THEN
         SELECT CASE(TRIM(opt_type))
            CASE('lmdif')
               WRITE(6,*) '    OPTIMIZER: Levenberg-Mardquardt'
               WRITE(6,*) '    NFUNC_MAX: ',nfunc_max
               WRITE(6,'(A,2X,1ES12.4)') '         FTOL: ',ftol
               WRITE(6,'(A,2X,1ES12.4)') '         XTOL: ',xtol
               WRITE(6,'(A,2X,1ES12.4)') '         GTOL: ',gtol
               WRITE(6,'(A,2X,1ES12.4)') '       EPSFCN: ',epsfcn
               WRITE(6,*) '         MODE: ',mode
               WRITE(6,*) '       FACTOR: ',factor
            CASE('lmdif_bounded')
               WRITE(6,*) '    OPTIMIZER: Levenberg-Mardquardt (Bounded)'
               WRITE(6,*) '    NFUNC_MAX: ',nfunc_max
               WRITE(6,'(A,2X,1ES12.4)') '         FTOL: ',ftol
               WRITE(6,'(A,2X,1ES12.4)') '         XTOL: ',xtol
               WRITE(6,'(A,2X,1ES12.4)') '         GTOL: ',gtol
               WRITE(6,'(A,2X,1ES12.4)') '       EPSFCN: ',epsfcn
               WRITE(6,*) '         MODE: ',mode
               WRITE(6,*) '       FACTOR: ',factor
            CASE('eval_xvec')
               WRITE(6,*) '    OPTIMIZER: XVEC Evlauation'
               WRITE(6,*) '    FILE:     ',TRIM(xvec_file)
            CASE('one_iter','single','eval','single_iter')
               WRITE(6,*) '    OPTIMIZER: SINGLE_ITERATION'
               WRITE(6,*) '    NFUNC_MAX: ',nfunc_max
            CASE('one_iter_norm')
               WRITE(6,*) '    OPTIMIZER: SINGLE_ITERATION FOR NORMALIZTION'
            CASE('gade')
               WRITE(6,*) '    OPTIMIZER: Differential Evolution'
               WRITE(6,*) '         NPOP: ',npopulation
               WRITE(6,*) '    NFUNC_MAX: ',nfunc_max
               WRITE(6,*) '       FACTOR: ',factor,'  (Mutation Scaling Factor)'
               WRITE(6,*) '       EPSFCN: ',epsfcn,'  (Crossover Factor)'
               WRITE(6,*) '         MODE: ',mode,'  (Strategy)'
               WRITE(6,*) '  CR_STRATEGY: ',cr_strategy
               IF(lrestart) WRITE(6,*) ' Restart file: ',TRIM('gade_restart.'//TRIM(id_string))
            CASE('map')
               WRITE(6,*) '    OPTIMIZER: Parameter Space Mapping'
               WRITE(6,*) '         NPOP: ',npopulation
               WRITE(6,*) '         NOPT: ',noptimizers
               WRITE(6,*) '         NDIV: ',mode
               WRITE(6,*) '            N: ',nvars
               WRITE(6,*) '            M: ',mtargets
               WRITE(6,*) '        NFUNC: ',mode**nvars
            CASE('map_linear')
               WRITE(6,*) '    OPTIMIZER: Linear Mapping'
               WRITE(6,*) '         NOPT: ',noptimizers
               WRITE(6,*) '         NDIV: ',mode
               WRITE(6,*) '            N: ',nvars
               WRITE(6,*) '            M: ',mtargets
               WRITE(6,*) '        NFUNC: ',mode*mode
            CASE('map_plane')
               WRITE(6,*) '    OPTIMIZER: Hyperplane Mapping'
               WRITE(6,*) '         NOPT: ',noptimizers
               WRITE(6,*) '         NDIV: ',mode
               WRITE(6,*) '       FACTOR: ',factor,'  (Scale factor for vectors)'
               WRITE(6,*) '            N: ',nvars
               WRITE(6,*) '            M: ',mtargets
               WRITE(6,*) '        NFUNC: ',mode*mode
            CASE('map_hypers')
               WRITE(6,*) '    OPTIMIZER: Hypersphere Mapping'
               WRITE(6,*) '         NRAD: ',nfunc_max
               WRITE(6,*) '         NPOP: ',npopulation
               WRITE(6,*) '         NOPT: ',noptimizers
               WRITE(6,*) '         DRHO: ',factor
               WRITE(6,*) '         ERHO: ',epsfcn
               WRITE(6,*) '            N: ',nvars
               WRITE(6,*) '            M: ',mtargets
            CASE('pso')
               WRITE(6,*) '    OPTIMIZER: Particle Swarm'
               WRITE(6,'(A,2X,1ES12.4)') '         FTOL: ',ftol
               WRITE(6,'(A,2X,1ES12.4)') '         XTOL: ',xtol
               WRITE(6,'(A,2X,1I5)')     '     NFUNC_MAX: ',nfunc_max
               WRITE(6,'(A,2X,1ES12.4)') '       C_local: ',epsfcn
               WRITE(6,'(A,2X,1ES12.4)') '      C_global: ',1.0
               WRITE(6,'(A,2X,1ES12.4)') '        Vscale: ',factor
               WRITE(6,'(A,2X,1I5)')     '          NPOP: ',npopulation
            CASE('rocket')
               WRITE(6,*) '    OPTIMIZER: Rocket'
               WRITE(6,'(A,2X,1ES12.4)') '         FTOL: ',ftol
               WRITE(6,'(A,2X,1ES12.4)') '         XTOL: ',xtol
               WRITE(6,'(A,2X,1I5)')     '     NFUNC_MAX: ',nfunc_max
               WRITE(6,'(A,2X,1ES12.4)') '       C_local: ',c1
               WRITE(6,'(A,2X,1ES12.4)') '      C_global: ',c2
               WRITE(6,'(A,2X,1ES12.4)') '        Vscale: ',factor
               WRITE(6,'(A,2X,1I5)')     '          NPOP: ',npopulation
         CASE DEFAULT
            ! This 'case default' section is entered if either a MANGO algorithm is selected, 
            ! or if an invalid opt_type is selected. In the latter case
            ! there is no need to print a warning here, since the next select case block
            ! will catch it.
         END SELECT
         IF (lauto_domain) WRITE(6,*) '  !!!!!! AUTO_DOMAIN Calculation !!!!!!!'
      END IF

      ! DEFAULT
      lskip_min = .false.
      ldeleteopt = .TRUE.
      lsocleanup = .TRUE.

      ! Do runs
      SELECT CASE(TRIM(opt_type))
         CASE('lmdif')
            ALLOCATE(ipvt(nvars))
            ALLOCATE(qtf(nvars),wa1(nvars),wa2(nvars),wa3(nvars),&
                     wa4(mtargets),fvec(mtargets))
            ALLOCATE(fjac(mtargets,nvars))
            fvec     = 0.0
            nprint   = 0
            info     = 0
            nfev     = 0
            ldfjac   = mtargets
            vars_min = -bigno; vars_max = bigno
            WHERE(vars > bigno) vars_max = 1E30
            CALL lmdif(stellopt_fcn, mtargets, nvars, vars, fvec, &
                       ftol, xtol, gtol, nfunc_max, epsfcn, diag, mode, &
                       factor, nprint, info, nfev, fjac, ldfjac, ipvt, &
                       qtf, wa1, wa2, wa3, wa4,vars_min,vars_max)
            DEALLOCATE(ipvt, qtf, wa1, wa2, wa3, wa4, fjac)
         CASE('lmdif_bounded')
            ALLOCATE(ipvt(nvars))
            ALLOCATE(qtf(nvars),wa1(nvars),wa2(nvars),wa3(nvars),&
                     wa4(mtargets),fvec(mtargets))
            ALLOCATE(fjac(mtargets,nvars))
            fvec     = 0.0
            nprint   = 0
            info     = 0
            nfev     = 0
            ldfjac   = mtargets
            CALL lmdif(stellopt_fcn, mtargets, nvars, vars, fvec, &
                       ftol, xtol, gtol, nfunc_max, epsfcn, diag, mode, &
                       factor, nprint, info, nfev, fjac, ldfjac, ipvt, &
                       qtf, wa1, wa2, wa3, wa4,vars_min,vars_max)
            DEALLOCATE(ipvt, qtf, wa1, wa2, wa3, wa4, fjac)
         CASE('eval_xvec')
            ldeleteopt = .FALSE.
            lskip_min = .true.
            CALL xvec_eval(stellopt_fcn,nvars,mtargets,xvec_file)
         CASE('one_iter','single','eval','single_iter')
            lskip_min = .true.
            ldeleteopt = .FALSE.
            ALLOCATE(fvec(mtargets))
            fvec     = 0.0
            info     = FLAG_SINGLETASK
            nfev     = 0
            CALL stellopt_fcn(mtargets, nvars, vars,fvec,info, nfev)
            c1 = enorm(mtargets,fvec)
            iunit = 12; info = 0
            CALL safe_open(iunit,info,'xvec.dat','unknown','formatted',ACCESS_IN='APPEND')
            WRITE(iunit,'(2(2X,I5.5))') nvars,0
            WRITE(iunit,'(10ES22.12E3)') vars(1:nvars)
            WRITE(iunit,'(ES22.12E3)') c1
            CLOSE(iunit)
            info = FLAG_CLEANUP
            IF (myid == master) info = flag_cleanup_lev
            call stellopt_fcn(mtargets, nvars, vars, fvec, info, nfev)
         CASE('one_iter_norm')
            ldeleteopt = .FALSE.
            lsocleanup = .FALSE.
            ALLOCATE(fvec(mtargets))
            fvec     = 0.0
            info     = FLAG_SINGLETASK
            nfev     = 0
            IF (myid == master) CALL stellopt_fcn(mtargets, nvars, vars,fvec,info, nfev)
         CASE('gade')
            ALLOCATE(fvec(mtargets))
         !   Notes on this
         !   factor: Mutation Scaling Factor
         !   epsfcn: Crossover Factor
         !   mode:   strategy 
         !   cr_strategy:  (0: exponential, 1: binomial)
            target_fitness = 0.0         ! Desired Fitness (chisq)
            npop           = npopulation ! Population Size (10*nvars is good)
            nprint         = 1           ! Keep every minimum
            iunit          = 27          ! Eventually we want to reinstate this with iunit=6
            iunit_restart  = 28
            IF (myid == master) THEN
               CALL safe_open(iunit,info,TRIM('gade_data.'//TRIM(id_string)),'unknown','formatted',ACCESS_IN='APPEND')
               INQUIRE(FILE=TRIM('gade_restart.'//TRIM(id_string)),EXIST=lfile_exists)
               IF (lfile_exists) THEN
                  CALL safe_open(iunit_restart,info,TRIM('gade_restart.'//TRIM(id_string)),'old','formatted')
               ELSE
                  CALL safe_open(iunit_restart,info,TRIM('gade_restart.'//TRIM(id_string)),'unknown','formatted')
               END IF
            END IF
            nfev           = 0
            CALL DE2_Evolve(stellopt_fcn,mtargets,nvars,npopulation,&
                            vars_min,vars_max,vars,fvec,nfunc_max,&
                            factor,epsfcn,mode,cr_strategy,iunit,&
                            iunit_restart,lrestart)
            CLOSE(iunit)
            CLOSE(iunit_restart)
         CASE('map')
            ldeleteopt = .FALSE.
            lskip_min = .true.
            nprint = 6
            lno_restart = .true.
            CALL MAP(stellopt_fcn,nvars,mtargets,vars_min,vars_max,npopulation,nprint,mode,MPI_COMM_STEL)
            info = 5
            IF (lverb) THEN
                WRITE(6,*) '!!!!!  PARAMETER SPACE MAPPING DONE  !!!!!'
                WRITE(6,*) '       See map.dat for data               '
            END IF
         CASE('map_linear')
            ldeleteopt = .FALSE.
            lskip_min = .true.
            nprint = 6
            lno_restart = .true.
            CALL MAP_LINEAR(stellopt_fcn,nvars,mtargets,vars,vars_min,vars_max,noptimizers,nprint,mode)
            info = 5
            IF (lverb) THEN
                WRITE(6,*) '!!!!!  LINEAR MAPPING DONE  !!!!!'
                WRITE(6,*) '       See map.dat for data               '
            END IF
         CASE('map_plane')
            ldeleteopt = .FALSE.
            lskip_min = .true.
            nprint = 6
            lno_restart = .true.
            CALL MAP_PLANE(stellopt_fcn,nvars,mtargets,vars,vars_min,vars_max,factor,noptimizers,nprint,mode)
            info = 5
            IF (lverb) THEN
                WRITE(6,*) '!!!!!  HYPERPLANE MAPPING DONE  !!!!!'
                WRITE(6,*) '       See map.dat for data               '
            END IF
         CASE('map_hypers')
            ldeleteopt = .FALSE.
            lskip_min = .true.
            ALLOCATE(wa1(nvars),fvec(mtargets))
            wa1 = vars
            nprint = 6
            lno_restart = .true.
            CALL MAP_HYPERS(stellopt_fcn,mtargets,nvars,noptimizers,vars_min,vars_max,&
                            wa1,fvec,factor,epsfcn,nfunc_max,MPI_COMM_STEL)
            IF (lverb) THEN
                WRITE(6,*) '!!!!!  HYPERSHPERE MAPPING DONE  !!!!!'
            END IF
            DEALLOCATE(wa1)
         CASE('pso')
            ALLOCATE(wa1(nvars),fvec(mtargets))
            wa1 = vars
            lno_restart = .TRUE.
            c1 = epsfcn
            c2 = 1.0
            CALL PSO_Evolve(stellopt_fcn,mtargets,nvars,npopulation,vars_min,vars_max,&
                            wa1,fvec,c1,c2,factor,ftol,xtol,nfunc_max)
            DEALLOCATE(wa1)
         CASE('rocket')
            ALLOCATE(wa1(nvars))
            wa1 = vars
            lno_restart = .TRUE.
            c1 = epsfcn
            c2 = 1.0
            CALL ROCKET_Evolve(stellopt_fcn,mtargets,nvars,npopulation,vars_min,vars_max,&
                            wa1,fvec,c1,c2,factor,ftol,xtol,nfunc_max)
            DEALLOCATE(wa1)
         CASE DEFAULT
            CALL stellopt_optimize_mango(used_mango_algorithm,nfev)
            IF (.not. used_mango_algorithm) THEN
               WRITE(6,*) '!!!!!  UNKNOWN OPTIMIZATION TYPE  !!!!!'
               WRITE(6,*) '       OPT_TYPE: ',TRIM(opt_type)
               WRITE(6,*)
               RETURN
            END IF
      END SELECT
      ! Now output the minimum files
      IF (myid == master .and. info /= 5 .and. .not.lskip_min) THEN
         IF (.NOT.ALLOCATED(fvec)) ALLOCATE(fvec(mtargets))
         IF (lrenorm) CALL stellopt_renorm(mtargets,fvec)
         nfev_save = nfev
         ier=0
         CALL stellopt_fcn(mtargets,nvars,vars,fvec,ier,nfev)
         nfev = nfev_save
         ier=-500
         CALL stellopt_fcn(mtargets,nvars,vars,fvec,ier,nfev)
      END IF
!DEC$ IF DEFINED (STELZIP)
      IF (myid == master .and. lsocleanup) THEN
         IF (ldeleteopt) THEN
            ! Remove the *_opt* files
            WRITE(6,*) ' Cleaning up _opt files'; CALL FLUSH(6); ier = 0; ierr_mpi = 0; cmdtxt=''
            CALL EXECUTE_COMMAND_LINE("rm -rf *_opt*",WAIT=.TRUE.,EXITSTAT=ier,CMDSTAT=ierr_mpi,CMDMSG=cmdtxt)
            WRITE(6,*) ' rm: EXITSTAT=',ier,' CMDSTAT=',ierr_mpi; CALL FLUSH(6)
            WRITE(6,*) '     MESSAGE: ',TRIM(cmdtxt); CALL FLUSH(6)
         END IF
         ! Zip up the results
         WRITE(6,*) ' Zipping files'; CALL FLUSH(6); ier = 0; ierr_mpi = 0; cmdtxt=''
         CALL EXECUTE_COMMAND_LINE("zip -r stellopt_files.zip *",WAIT=.TRUE.,EXITSTAT=ier,CMDSTAT=ierr_mpi,CMDMSG=cmdtxt)
         WRITE(6,*) ' zip: EXITSTAT=',ier,' CMDSTAT=',ierr_mpi; CALL FLUSH(6)
         WRITE(6,*) '     MESSAGE: ',TRIM(cmdtxt); CALL FLUSH(6)
         ier = 0; ierr_mpi=0
      END IF
!DEC$ ENDIF
      IF (ALLOCATED(fvec)) DEALLOCATE(fvec)

!DEC$ IF DEFINED (MPI_OPT)
      CALL MPI_BARRIER(MPI_COMM_STEL,ierr_mpi)                       !PPPL -SAL
      IF (ierr_mpi /= MPI_SUCCESS) CALL handle_err(MPI_ERR,'stellopt_optimize2',ierr_mpi)
!DEC$ ENDIF
        
!----------------------------------------------------------------------
!     END SUBROUTINE
!----------------------------------------------------------------------
      END SUBROUTINE stellopt_optimize
