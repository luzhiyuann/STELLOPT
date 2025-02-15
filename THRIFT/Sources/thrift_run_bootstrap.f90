!-----------------------------------------------------------------------
!     Subroutine:    thrift_run_bootstrap
!     Authors:       L. van Ham
!     Date:          11/XX/2022
!     Description:   This subroutine calls the relevant bootstrap
!                    current code.
!-----------------------------------------------------------------------
SUBROUTINE thrift_run_bootstrap
   !-----------------------------------------------------------------------
   !     Libraries
   !-----------------------------------------------------------------------
   USE thrift_runtime
   USE thrift_vars
   USE thrift_equil
   USE thrift_profiles_mod
   USE thrift_funcs
   USE booz_params, ONLY: lsurf_boz
   !-----------------------------------------------------------------------
   !     Local Variables
   !        ier         Error flag
   !-----------------------------------------------------------------------
   IMPLICIT NONE
   INTEGER :: i, ier
   REAL(rprec) :: epsilon, pp1, pm1, ds, rho, s
   REAL(rprec), DIMENSION(:), ALLOCATABLE ::  j_temp
   !----------------------------------------------------------------------
   !     BEGIN SUBROUTINE
   !----------------------------------------------------------------------

   ier = 0

   ! Just always check it's deallocated      
   IF (ALLOCATED(lsurf_boz)) DEALLOCATE(lsurf_boz)

   ! Check to make sure we're not zero beta
   IF (eq_beta == 0) THEN
      THRIFT_JBOOT(:,mytimestep) = 0
      RETURN
   END IF

   SELECT CASE(TRIM(bootstrap_type))
      CASE ('model','simple','test')

         ! j_BS = sqrt(epsilon) Rmajor *dp/dPhi
         ! epsilon = a/R (inverse aspect ratio)
         ! dp/dPhi : Pa/Wb (toroidal flux derivative)
         ! dp/dPhi = dp/ds * ds/dPhi
         !         = dp/ds / Phi_edge
         ! j_BS = s*sqrt(epsilon)*Rmajor/Phi_edge*dp/ds
         !      = s*sqrt(aminor*Rmajor)/Phi_edge*dp/ds
         
         epsilon = eq_Aminor/eq_Rmajor ! Inverse Aspect Ratio
         THRIFT_JBOOT(:,mytimestep) = SQRT(eq_Aminor*eq_Rmajor)/eq_phiedge*THRIFT_S*THRIFT_PPRIME(:,mytimestep)

      CASE ('bootsj')
         ALLOCATE(lsurf_boz(ns_eq))
         lsurf_boz = .FALSE.
         lsurf_boz(2:ns_eq) = .TRUE.
         CALL thrift_paraexe('booz_xform',proc_string,lscreen_subcodes)
         CALL thrift_paraexe('bootsj',proc_string,lscreen_subcodes)
      CASE('read_from_file')
         DO i = 1, nsj
            s = THRIFT_S(i)
            rho = SQRT(s)
            CALL get_prof_JBS(rho,THRIFT_T(mytimestep),THRIFT_JBOOT(i,mytimestep))
         END DO
      CASE ('dkespenta')
         ALLOCATE(lsurf_boz(ns_eq))
         lsurf_boz = .FALSE.
         DO i = 1, DKES_NS_MAX
            IF (DKES_K(i)<1) CYCLE
            lsurf_boz(DKES_K(i)) = .TRUE.
         END DO
         CALL thrift_paraexe('booz_xform',proc_string,lscreen_subcodes)
         CALL thrift_paraexe('dkes',proc_string,lscreen_subcodes)
         CALL thrift_paraexe('penta',proc_string,lscreen_subcodes)
      CASE ('sfincs')
   END SELECT

   RETURN
   !----------------------------------------------------------------------
   !     END SUBROUTINE
   !----------------------------------------------------------------------
   END SUBROUTINE thrift_run_bootstrap

