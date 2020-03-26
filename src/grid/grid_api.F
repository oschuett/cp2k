!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
MODULE grid_api
   USE grid_collocate,                  ONLY: &
        GRID_FUNC_AB, GRID_FUNC_ADBmDAB_X, GRID_FUNC_ADBmDAB_Y, GRID_FUNC_ADBmDAB_Z, &
        GRID_FUNC_ARDBmDARB_XX, GRID_FUNC_ARDBmDARB_XY, GRID_FUNC_ARDBmDARB_XZ, &
        GRID_FUNC_ARDBmDARB_YX, GRID_FUNC_ARDBmDARB_YY, GRID_FUNC_ARDBmDARB_YZ, &
        GRID_FUNC_ARDBmDARB_ZX, GRID_FUNC_ARDBmDARB_ZY, GRID_FUNC_ARDBmDARB_ZZ, &
        GRID_FUNC_DABpADB_X, GRID_FUNC_DABpADB_Y, GRID_FUNC_DABpADB_Z, GRID_FUNC_DADB, &
        GRID_FUNC_DX, GRID_FUNC_DXDX, GRID_FUNC_DXDY, GRID_FUNC_DY, GRID_FUNC_DYDY, &
        GRID_FUNC_DYDZ, GRID_FUNC_DZ, GRID_FUNC_DZDX, GRID_FUNC_DZDZ, collocate_pgf_product
   USE grid_integrate,                  ONLY: integrate_pgf_product
   USE message_passing,                 ONLY: mp_sum
   USE ISO_C_BINDING,                   ONLY: C_INT, C_FUNPTR, C_FUNLOC, C_NULL_CHAR, C_CHAR, C_LONG
   USE OMP_LIB, ONLY: omp_get_thread_num
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'grid_api'

   PUBLIC :: GRID_FUNC_AB
   PUBLIC :: GRID_FUNC_DADB
   PUBLIC :: GRID_FUNC_ADBmDAB_X
   PUBLIC :: GRID_FUNC_ADBmDAB_Y
   PUBLIC :: GRID_FUNC_ADBmDAB_Z
   PUBLIC :: GRID_FUNC_ARDBmDARB_XX
   PUBLIC :: GRID_FUNC_ARDBmDARB_XY
   PUBLIC :: GRID_FUNC_ARDBmDARB_XZ
   PUBLIC :: GRID_FUNC_ARDBmDARB_YX
   PUBLIC :: GRID_FUNC_ARDBmDARB_YY
   PUBLIC :: GRID_FUNC_ARDBmDARB_YZ
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZX
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZY
   PUBLIC :: GRID_FUNC_ARDBmDARB_ZZ
   PUBLIC :: GRID_FUNC_DABpADB_X
   PUBLIC :: GRID_FUNC_DABpADB_Y
   PUBLIC :: GRID_FUNC_DABpADB_Z
   PUBLIC :: GRID_FUNC_DX
   PUBLIC :: GRID_FUNC_DY
   PUBLIC :: GRID_FUNC_DZ
   PUBLIC :: GRID_FUNC_DXDY
   PUBLIC :: GRID_FUNC_DYDZ
   PUBLIC :: GRID_FUNC_DZDX
   PUBLIC :: GRID_FUNC_DXDX
   PUBLIC :: GRID_FUNC_DYDY
   PUBLIC :: GRID_FUNC_DZDZ

   PUBLIC :: collocate_pgf_product
   PUBLIC :: integrate_pgf_product
   PUBLIC :: grid_library_init
   PUBLIC :: grid_library_finalize

CONTAINS

! **************************************************************************************************
!> \brief Initialize grid library
! **************************************************************************************************
   SUBROUTINE grid_library_init()
      INTERFACE
         SUBROUTINE grid_globals_init(thread_num) BIND(C, name="grid_globals_init")
            IMPORT :: C_INT
            INTEGER(KIND=C_INT), VALUE                :: thread_num
         END SUBROUTINE grid_globals_init
      END INTERFACE

      ! Passing thread_num to ensure that Fortran and C runtimes agree on OpenMP thread numbers.

!$OMP PARALLEL DEFAULT(NONE)
      CALL grid_globals_init(thread_num=omp_get_thread_num())
!$OMP END PARALLEL

   END SUBROUTINE grid_library_init

! **************************************************************************************************
!> \brief Finalize grid library
! **************************************************************************************************
   SUBROUTINE grid_library_finalize(mpi_comm, unit_nr)
     INTEGER, INTENT(IN)                             :: mpi_comm, unit_nr

     INTERFACE
        SUBROUTINE grid_globals_finalize(mpi_sum_func, print_func) BIND(C, name="grid_globals_finalize")
           IMPORT :: C_FUNPTR
           TYPE(C_FUNPTR), VALUE                     :: mpi_sum_func
           TYPE(C_FUNPTR), VALUE                     :: print_func
        END SUBROUTINE grid_globals_finalize
     END INTERFACE

     ! Since Fortran unit and mpi-group numbers can't be used from C, we pass function pointers instead.

!$OMP PARALLEL DEFAULT(NONE)
     CALL grid_globals_finalize(mpi_sum_func=C_FUNLOC(mpi_sum_func), &
                                print_func=C_FUNLOC(print_func))
!$OMP END PARALLEL

     CONTAINS

     SUBROUTINE mpi_sum_func(number) BIND(C)
        INTEGER(KIND=C_LONG), INTENT(INOUT)                :: number
        CALL mp_sum(number, mpi_comm)
     END SUBROUTINE mpi_sum_func

     SUBROUTINE print_func(message) BIND(C)
        CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN)          :: message(*)
        CHARACTER(LEN=1000)                                :: buffer
        INTEGER                                            :: i
        IF (unit_nr > 0) THEN
           ! Convert C char array into Fortran string.
           buffer = ""
           DO i = 1, LEN(buffer)
              IF (message(i) == C_NULL_CHAR) EXIT
              buffer(i:i) = message(i)
           END DO
           ! Print the message.
           WRITE (unit_nr, FMT="(A)", ADVANCE="NO") buffer(1:i-1)
        END IF
      END SUBROUTINE print_func
   END SUBROUTINE grid_library_finalize

END MODULE grid_api
