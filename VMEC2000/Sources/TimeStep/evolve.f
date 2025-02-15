      SUBROUTINE evolve(time_step, ier_flag, liter_flag, lscreen)
      USE vmec_main
      USE vmec_params, ONLY: bad_jacobian_flag, successful_term_flag,
     &                       norm_term_flag, time_control_flag
      USE xstuff
      USE precon2d, ONLY: ictrl_prec2d, l_comp_prec2D, 
     &                    compute_blocks_par, compute_blocks
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: ZeroLastNType, CopyLastNtype, 
     &                                SaxpbyLastNtype, CompareEdgeValues
      USE timer_sub
      USE vmec_params, ONLY: ntmax
      USE gmres_mod
!  Comment Out below JDH 2010-08-03
!      USE vmec_history
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp)            :: time_step          !, r0dot
      INTEGER, INTENT(INOUT) :: ier_flag
      LOGICAL, INTENT(INOUT) :: liter_flag
      LOGICAL, INTENT(IN)    :: lscreen
!-----------------------------------------------``
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CHARACTER(LEN=*), PARAMETER :: fcn_message =
     &   'External calls to FUNCT3D: '
!      REAL(dp), PARAMETER :: r0dot_threshold = 5.E-06_dp
      REAL(dp) :: fsq1, dtau, b1, bprec, fac
      LOGICAL :: lfinal_mesh
      INTEGER :: lcount
      INTEGER, SAVE :: iter_on
      REAL(dp) :: f3dt1, f3dt2, tevon, tevoff

C-----------------------------------------------
!     IF TROUBLE CONVERGING, TRY TO RECOMPUTE PRECONDITIONER ONCE MORE...
!      IF (ictrl_prec2d.eq.1 .and. iter2.eq.(iter_on+40)) 
!     1     ictrl_prec2d = 0

!  JDH 2011-09-15 Add condition to lfinal_mesh, that iter2 - iter1 > 5
!    (5 was picked out of a hat)
!    Purpose is to keep preconditioning from being turned on immediately upon
!    a restart (V3FIT)
!   The final .and. clause is a bit complicated. The purpose of the
!   .not. lv3fit is so that if v3fita is not running, the final .and. clause
!   is always true. Read as "and, if lv3fit, then must also have iter2 - iter1 > 5"
!      lfinal_mesh = (ns .eq. ns_maxval) .and. (ictrl_prec2d.eq.0)
!     1              .and. (itype_precon.ne.0)

      CALL second0(tevon)


      lfinal_mesh = ns           .EQ. ns_maxval .and.
     &              ictrl_prec2d .EQ. 0         .and.
     &              itype_precon .ne. 0         .and.
     &              (.not.l_v3fit .or. iter2 - iter1 .ge. 5)

      IF (iter2 .lt. 10) THEN
         ictrl_prec2d = 0
         lqmr = .false.
         iter_on = -1
      ELSE IF (lfinal_mesh .and. 
     &         fsqr + fsqz + fsql .lt. prec2d_threshold) THEN
         lqmr = (itype_precon .GE. 2)
         lfirst = (lqmr .AND. iter_on.EQ.-1) 

!
!        INITIATES 2D PRECONDITIONER CALCULATION
!
         IF (iter_on .EQ. -1) THEN
            IF (pre_niter .eq. -1) THEN
               IF (lqmr) THEN
                  nstep = 5
                  niter = iter2+100                   !Limit # preconditioner steps
               ELSE
                  nstep = 20
                  niter = iter2+400
               END IF

               IF (rank .eq. 0) THEN
                  WRITE (6,1000) niter
                  WRITE (nthreed,1000) niter
               END IF
            ELSE IF (pre_niter .ge. 0) THEN
               niter = iter2 + pre_niter

               IF (rank .eq. 0) THEN
                  WRITE (6,1000) niter
                  WRITE (nthreed,1000) niter
               END IF
            END IF

            iter_on = iter2                        !Flag to monitor progress of preconditioner
         ELSE
            iter_on = iter2-11
         END IF

!SPH022111: ADD NEW CONTROL PARAMETER, l_comp_prec2D, TO FORCE RECALCULATION
!           OF PRECONDITIONING BLOCKS IN V3FIT, FOR EXAMPLE
         IF (lfirst .OR. l_comp_prec2D) THEN
            IF (l_v3fit) THEN
               WRITE(*,*) 'VMEC Evolve:compute_blocks'
            END IF
            IF (PARVMEC) THEN
               CALL compute_blocks_par (pxc,pxcdot,pgc)
            ELSE
               CALL compute_blocks (xc,xcdot,gc)
            END IF
         END IF
         IF(l_v3fit) WRITE(*,*) 'VMEC Evolve:prec2d_On iter2 =', iter2
         l_comp_prec2D = .FALSE.
         ictrl_prec2d = 1
         time_step = 0.50_dp
         iter1 = iter2-1; fsq = fsqr1 + fsqz1 + fsql1

         IF (PARVMEC) THEN
            CALL CopyLastNtype(pxstore, pxc)
            CALL ZeroLastNType(pxcdot)
         ELSE
            xc = xstore
            xcdot = 0
         END IF
      END IF
!
!     CHECK NS_RESLTN
!
      IF (NS_RESLTN > num_grids) THEN
         ier_flag = time_control_flag
         RETURN
      ENDIF

!
!     COMPUTE MHD FORCES
!     MUST CALL funct3d EVEN WHEN IN 2D PRECONDITIONING MODE, SINCE
!     INITIAL RESIDUALS MUST BE KNOWN WHEN CALLING gmres_fun, etc.
!
      CALL second0(f3dt1)
      f3d_num(NS_RESLTN) = f3d_num(NS_RESLTN)+1
      IF (PARVMEC) THEN
         CALL funct3d_par(lscreen, ier_flag)
      ELSE
         CALL funct3d(lscreen, ier_flag)
      END IF
      CALL second0(f3dt2)
      f3d_time(NS_RESLTN) = f3d_time(NS_RESLTN) + (f3dt2 - f3dt1)
      funct3d_time = funct3d_time + (f3dt2 - f3dt1)

!
!     COMPUTE ABSOLUTE STOPPING CRITERION
      IF (iter2.EQ.1 .and. irst.EQ.2) THEN
         ier_flag = bad_jacobian_flag
         RETURN
!  JDH 2012-04-24. Revise this absolute stopping criterion, so that if v3fit
!    is running, then have to iterate at least 2 * nvacskip steps
!    (2 picked out of a hat) (nvacskip - to make sure vacuum gets updated)
!    before returning.

      ELSE IF (fsqr .le. ftolv .and.
     &         fsqz .le. ftolv .and.
     &         fsql .le. ftolv) THEN
         liter_flag = .false.
         ier_flag = successful_term_flag
         RETURN
      ENDIF


!SPH:042117: MOVE TIME STEP CONTROL HERE (FROM END OF EQSOLVE) TO AVOID
!STORING A POSSIBLE irst=2 STATE
      CALL TimeStepControl(ier_flag, PARVMEC)

      IF (lqmr) THEN
         IF (PARVMEC) THEN
            CALL gmres_fun_par(ier_flag, itype_precon - 1, lscreen)
            IF (.NOT.lfreeb) CALL CompareEdgeValues(pxc, pxsave)
         ELSE
            CALL gmres_fun(ier_flag, itype_precon-1)
            IF (.NOT.lfreeb) THEN
               DO lcount = ns, 2*irzloff, ns
                  IF (xsave(lcount) .NE. xc(lcount)) THEN
                     PRINT *, ' xsave = ',xsave(lcount),' != xc = ',
     &                        xc(lcount),' for lcount = ',lcount
                  END IF
               END DO
            END IF
         END IF

         RETURN
      END IF

!     COMPUTE DAMPING PARAMETER (DTAU) AND EVOLVE
!     R, Z, AND LAMBDA ARRAYS IN FOURIER SPACE

      fsq1 = fsqr1 + fsqz1 + fsql1

      IF (iter2 .EQ. iter1) otau(:ndamp) = cp15/time_step

      IF (ictrl_prec2d .EQ. 0) THEN
         bprec = 1
      ELSE
         bprec = 6
      END IF

      dtau = bprec*cp15
      IF (iter2 .GT. iter1 .AND.
     &    fsq1*fsq .NE. zero) THEN
         dtau = MIN(ABS(LOG(fsq1/fsq)), dtau)
      END IF

      fsq = fsq1

      otau(1:ndamp-1) = otau(2:ndamp)

      IF (iter2 .GT. iter1) otau(ndamp) = dtau/time_step
!REMOVED 071505: OTHERWISE I=1 STATE REPEATED (SKIP THIS TO GET OUT OF ITER2=1 STATE)
!     IF (iter2 .le. 1) RETURN

      otav = SUM(otau(:ndamp))/ndamp
      dtau = time_step*otav/2

      b1  = one - dtau
      fac = one/(one + dtau)

!
!     THIS IS THE TIME-STEP ALGORITHM. IT IS ESSENTIALLY A CONJUGATE
!     GRADIENT METHOD, WITHOUT THE LINE SEARCHES (FLETCHER-REEVES),
!     BASED ON A METHOD GIVEN BY P. GARABEDIAN

!

      IF (PARVMEC) THEN
         IF (lactive) THEN
            CALL SaxpbyLastNtype(fac*time_step, pgc, fac*b1, pxcdot,
     &                           pxcdot)
            CALL SaxpbyLastNtype(time_step, pxcdot, one, pxc, pxc)
         END IF
      ELSE
         xcdot = fac*(b1*xcdot + time_step*gc)
         xc    = xc + time_step*xcdot
      END IF

      CALL second0(tevoff)
      evolve_time = evolve_time + (tevoff - tevon)

1000  FORMAT(2x,'Resetting the number of niter to ',i6)

      END SUBROUTINE evolve


      SUBROUTINE TimeStepControl(ier_flag, PARVMEC)
      USE vmec_main, ONLY: res0, res1, fsq, fsqr, fsqz, fsql,
     &                     irst, iter1, iter2, delt0r, dp
      USE vmec_params, ONLY: ns4, time_control_flag
      USE vparams, ONLY: c1pm2
      USE vmec_input, ONLY: nstep
      USE precon2d, ONLY: ictrl_prec2d
      USE parallel_include_module, ONLY: rank
      USE realspace
      IMPLICIT NONE
!     
!     STORES OR RETRIEVES XC STATE BASED ON IRST VALUE
!
      REAL(dp), PARAMETER :: fact = 1.E4_dp
      REAL(dp) :: fsq0
      INTEGER  :: ier_flag
      LOGICAL, INTENT(IN) :: PARVMEC

      ier_flag = 0
      fsq0 = fsqr+fsqz+fsql
      IF (iter2.EQ.iter1 .OR. res0.EQ.-1) THEN
         res0 = fsq
         res1 = fsq0
         CALL restart_iter(delt0r)
      END IF

      res0 = MIN(res0,fsq)
      res1 = MIN(res1,fsq0)


! Store current state (irst=1)
      IF (fsq.LE.res0 .AND. fsq0.LE.res1 .AND. irst.EQ.1) THEN 
         CALL restart_iter(delt0r)

      ELSE IF (ictrl_prec2d .NE. 0) THEN
         CALL restart_iter(delt0r)
         RETURN

      ELSE IF ((iter2-iter1) .GT. 10) THEN

! Residuals are growing in time, reduce time step
         IF (fsq.GT.fact*res0 .OR. fsq0.GT.fact*res1) THEN
            irst = 3
         END IF
      END IF

!     Retrieve previous good state
      IF (irst .NE. 1) THEN
         CALL restart_iter(delt0r)
         iter1 = iter2
         IF (PARVMEC) THEN
            CALL funct3d_par(.FALSE., ier_flag)
         ELSE
            CALL funct3d(.FALSE., ier_flag)
         END IF
         IF (irst .NE. 1 .and. irst .NE. 4) THEN
            ier_flag = time_control_flag
            !STOP 'Logic error in TimeStepControl!'
         END IF
      END IF

      END SUBROUTINE TimeStepControl
