!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Direct parametrization
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_direct
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_create, dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, &
        dbcsr_iterator_next_block, dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_set, dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE pao_potentials,                  ONLY: pao_guess_initial_potential
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE mathlib,                         ONLY: diamat_all
   USE pao_lnv,                         ONLY: pao_calc_grad_lnv_wrt_AB

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_direct'

   PUBLIC :: pao_param_init_direct, pao_param_finalize_direct, pao_calc_AB_direct
   PUBLIC :: pao_param_count_direct, pao_param_initguess_direct

CONTAINS

! **************************************************************************************************
!> \brief Initialize direct parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_direct(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      IF (pao%precondition) &
         CPABORT("PAO preconditioning not supported for selected parametrization.")

   END SUBROUTINE pao_param_init_direct

! **************************************************************************************************
!> \brief Finalize direct parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_direct()

      ! Nothing todo.

   END SUBROUTINE pao_param_finalize_direct

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_direct(qs_env, ikind, nparams)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_param_count_direct'

      INTEGER                                            :: pao_basis_size, pri_basis_size
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind(qs_kind_set(ikind), &
                       basis_set=basis_set, &
                       pao_basis_size=pao_basis_size)
      pri_basis_size = basis_set%nsgf

      nparams = pao_basis_size * pri_basis_size

   END SUBROUTINE pao_param_count_direct

! **************************************************************************************************
!> \brief Fills matrix_X with an initial guess
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_direct(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_param_initguess_direct'

      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X, my_block 
      INTEGER                                            :: handle, iatom, arow, acol, m, n, i
      TYPE(dbcsr_iterator_type)                          :: iter
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri

      CALL timeset(routineN, handle)

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,blk_sizes_pri,blk_sizes_pao) &
!$OMP PRIVATE(iter,arow,acol,iatom,n,m,block_X)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)

         n = blk_sizes_pri(iatom) ! size of primary basis
         m = blk_sizes_pao(iatom) ! size of pao basis
         ALLOCATE(my_block(n,m)) !TODO: do this without allocation
         my_block = 0.0_dp
         DO i = 1, m
            my_block(i,i) = 1.0_dp
         ENDDO
         block_X = my_block
         DEALLOCATE(my_block)

      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL timestop(handle)

   END SUBROUTINE pao_param_initguess_direct

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_direct(pao, qs_env, ls_scf_env, gradient)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_AB_direct'

      INTEGER                                            :: handle, acol, arow, iatom, n, m, i, j, k
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_X, block_A, block_B, AA, &
                                           AA_evecs, AA_inv, block_Ma, block_Mb, block_G
      REAL(dp), DIMENSION(:), POINTER                    :: AA_evals
      REAL(dp)                                           :: w
      LOGICAL                                            :: found
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_type)                                   :: matrix_Ma, matrix_Mb


      CALL timeset(routineN, handle)
      ls_mstruct => ls_scf_env%ls_mstruct

      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)
      ENDIF

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,matrix_Ma,matrix_Mb) &
!$OMP PRIVATE(iter,arow,acol,iatom,found,n,m,w,i,j,k) &
!$OMP PRIVATE(block_X,block_A,block_B,AA,AA_evals,AA_evecs,AA_inv) &
!$OMP PRIVATE(block_Ma, block_Mb, block_G)
      CALL dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_A, row=iatom, col=iatom, block=block_A, found=found)
         CPASSERT(ASSOCIATED(block_A))
         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_B, row=iatom, col=iatom, block=block_B, found=found)
         CPASSERT(ASSOCIATED(block_B))

         block_A = block_X
         n = SIZE(block_A, 1) ! size of primary basis
         m = SIZE(block_A, 2) ! size of pao basis

         WRITE(*,*) "OLE iatom: ", iatom, " block_X: ", block_X
         !WRITE(*,*) "OLE N: ", N, " M: ", M
         !DO i=1, n
         !  WRITE (*,*) "OLE A:", block_A(i,:)
         !ENDDO

         ALLOCATE(AA(m,m))
         AA = MATMUL(TRANSPOSE(block_A), block_A)

         !WRITE (*,*) "OLE AA:", AA

         ! diagonalize AA
         ALLOCATE (AA_evecs(m,m), AA_evals(m))
         AA_evecs(:,:) = AA
         CALL diamat_all(AA_evecs, AA_evals)

         !WRITE(*,*) "OLE evals: ", AA_evals

         ! build AA_inv
         ALLOCATE (AA_inv(m,m))
         AA_inv(:,:) = 0.0_dp
         DO k = 1, m
            w = 1.0_dp/SQRT(AA_evals(k))
            DO i = 1, m
            DO j = 1, m
               AA_inv(i, j) = AA_inv(i, j) + w * AA_evecs(i, k) * AA_evecs(j, k)
            ENDDO
            ENDDO
         ENDDO

         block_B = MATMUL(block_A, AA_inv)

         ! TURNING POINT (if calc grad) ------------------------------------------
         IF (gradient) THEN
            CALL dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
            CPASSERT(ASSOCIATED(block_G))
            CALL dbcsr_get_block_p(matrix=matrix_Ma, row=iatom, col=iatom, block=block_Ma, found=found)
            CALL dbcsr_get_block_p(matrix=matrix_Mb, row=iatom, col=iatom, block=block_Mb, found=found)
            ! don't check ASSOCIATED(block_M), it might have been filtered out.

            block_G(:,:) = 0.0_dp  !TODO: really needed?
            IF (ASSOCIATED(block_Ma)) THEN
               !WRITE(*,*) "OLE iatom: ", iatom, " block_Ma: ", block_Ma
               block_G = block_Ma
            ENDIF

            IF (ASSOCIATED(block_Mb)) THEN
               !WRITE(*,*) "OLE iatom: ", iatom, " block_Mb: ", block_Mb
               block_G = block_G + MATMUL(block_Mb, AA_inv) !TODO: just a wild guess and missing terms
            ENDIF

            !WRITE(*,*) "OLE iatom: ", iatom, " block_G: ", block_G

         ENDIF

         DEALLOCATE(AA, AA_evecs, AA_evals, AA_inv)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      IF (gradient) THEN
         CALL dbcsr_release(matrix_Ma)
         CALL dbcsr_release(matrix_Mb)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pao_calc_AB_direct

!! **************************************************************************************************
!!> \brief Calculate new matrix U and optionally its gradient G
!!> \param pao ...
!!> \param matrix_U ...
!!> \param matrix_M ...
!!> \param matrix_G ...
!! **************************************************************************************************
!   SUBROUTINE pao_calc_U_direct(pao, matrix_U, matrix_M, matrix_G)
!      TYPE(pao_env_type), POINTER                        :: pao
!      TYPE(dbcsr_type)                                   :: matrix_U
!      TYPE(dbcsr_type), OPTIONAL                         :: matrix_M, matrix_G
!
!      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_U_direct', routineP = moduleN//':'//routineN
!
!      COMPLEX(dp)                                        :: denom
!      COMPLEX(dp), DIMENSION(:), POINTER                 :: evalsblock_U,block_U0,block_X_full,evals,evecs
!      COMPLEX(dp), DIMENSION(:, :), POINTER              :: block_D, evecs
!      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, M, &
!                                                            N, nparams
!      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
!      LOGICAL                                            :: found
!      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_G_full, block_M, &
!                                                            block_tmp, block_U, block_U0, block_X, &
!                                                            block_X_full
!      TYPE(dbcsr_iterator_type)                          :: iter
!
!      CALL timeset(routineN, handle)
!
!      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)
!
!!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,matrix_M,matrix_G,blk_sizes_pri,blk_sizes_pao) &
!!$OMP PRIVATE(iter,arow,acol,iatom,N,M,nparams,i,j,k,found) &
!!$OMP PRIVATE(block_X,block_U,block_U0,block_X_full,evals,evecs) &
!!$OMP PRIVATE(block_M,block_G,block_D,block_tmp,block_G_full,denom)
!      CALL dbcsr_iterator_start(iter, pao%matrix_X)
!      DO WHILE (dbcsr_iterator_blocks_left(iter))
!         CALL dbcsr_iterator_next_block(iter, arow, acol, block_X)
!         iatom = arow; CPASSERT(arow == acol)
!         CALL dbcsr_get_block_p(matrix=matrix_U, row=iatom, col=iatom, block=block_U, found=found)
!         CPASSERT(ASSOCIATED(block_U))
!         CALL dbcsr_get_block_p(matrix=pao%matrix_U0, row=iatom, col=iatom, block=block_U0, found=found)
!         CPASSERT(ASSOCIATED(block_U0))
!
!         N = blk_sizes_pri(iatom) ! size of primary basis
!         M = blk_sizes_pao(iatom) ! size of pao basis
!         nparams = SIZE(block_X, 1)
!
!         ! block_X stores only rotations between occupied and virtuals
!         ! hence, we first have to build the full anti-symmetric exponent block
!         ALLOCATE (block_X_full(N, N))
!         block_X_full(:, :) = 0.0_dp
!         DO i = 1, nparams
!            block_X_full(MOD(i-1, M)+1, M+(i-1)/M+1) = +block_X(i, 1)
!            block_X_full(M+(i-1)/M+1, MOD(i-1, M)+1) = -block_X(i, 1)
!         ENDDO
!
!         ! diagonalize block_X_full
!         ALLOCATE (evals(N), evecs(N, N))
!         CALL diag_antisym(block_X_full, evecs, evals)
!
!         ! construct rotation matrix
!         block_U(:, :) = 0.0_dp
!         DO k = 1, N
!            DO i = 1, N
!               DO j = 1, N
!                  block_U(i, j) = block_U(i, j)+REAL(EXP(evals(k))*evecs(i, k)*CONJG(evecs(j, k)), dp)
!               ENDDO
!            ENDDO
!         ENDDO
!
!         block_U = MATMUL(block_U0, block_U) ! prepend initial guess rotation
!
!         ! TURNING POINT (if calc grad) ------------------------------------------
!         IF (PRESENT(matrix_G)) THEN
!            CPASSERT(PRESENT(matrix_M))
!
!            CALL dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
!            CPASSERT(ASSOCIATED(block_G))
!            CALL dbcsr_get_block_p(matrix=matrix_M, row=iatom, col=iatom, block=block_M, found=found)
!            ! don't check ASSOCIATED(block_M), it might have been filtered out.
!
!            ALLOCATE (block_D(N, N), block_tmp(N, N), block_G_full(N, N))
!            DO i = 1, N
!               DO j = 1, N
!                  denom = evals(i)-evals(j)
!                  IF (i == j) THEN
!                     block_D(i, i) = EXP(evals(i)) ! diagonal elements
!                  ELSE IF (ABS(denom) > 1e-10_dp) THEN
!                     block_D(i, j) = (EXP(evals(i))-EXP(evals(j)))/denom
!                  ELSE
!                     block_D(i, j) = 1.0_dp ! limit according to L'Hospital's rule
!                  ENDIF
!               ENDDO
!            ENDDO
!
!            IF (ASSOCIATED(block_M)) THEN
!               block_tmp = MATMUL(TRANSPOSE(block_U0), block_M)
!            ELSE
!               block_tmp = 0.0_dp
!            ENDIF
!            block_G_full = fold_derivatives(block_tmp, block_D, evecs)
!
!            ! return only gradient for rotations between occupied and virtuals
!            DO i = 1, nparams
!               block_G(i, 1) = 2.0_dp*block_G_full(MOD(i-1, M)+1, M+(i-1)/M+1)
!            ENDDO
!
!            DEALLOCATE (block_D, block_tmp, block_G_full)
!         ENDIF
!
!         DEALLOCATE (block_X_full, evals, evecs)
!
!      END DO
!      CALL dbcsr_iterator_stop(iter)
!!$OMP END PARALLEL
!
!      CALL timestop(handle)
!   END SUBROUTINE pao_calc_U_direct
!
! **************************************************************************************************
!> \brief Helper routine, for calculating derivatives
!> \param M ...
!> \param D ...
!> \param R ...
!> \return ...
! **************************************************************************************************
   FUNCTION fold_derivatives(M, D, R) RESULT(G)
      REAL(dp), DIMENSION(:, :), INTENT(IN)              :: M
      COMPLEX(dp), DIMENSION(:, :), INTENT(IN)           :: D, R
      REAL(dp), DIMENSION(SIZE(M, 1), SIZE(M, 1))        :: G

      CHARACTER(len=*), PARAMETER :: routineN = 'fold_derivatives', &
         routineP = moduleN//':'//routineN

      COMPLEX(dp), DIMENSION(:, :), POINTER              :: F, RF, RM, RMR
      INTEGER                                            :: n
      REAL(dp), DIMENSION(:, :), POINTER                 :: RFR

      n = SIZE(M, 1)

      ALLOCATE (RM(n, n), RMR(n, n), F(n, n), RF(n, n), RFR(n, n))

      RM = MATMUL(TRANSPOSE(CONJG(R)), TRANSPOSE(M))
      RMR = MATMUL(RM, R)
      F = RMR*D !Hadamard product
      RF = MATMUL(R, F)
      RFR = REAL(MATMUL(RF, TRANSPOSE(CONJG(R))))

      ! gradient dE/dX has to be anti-symmetric
      G = 0.5_dp*(TRANSPOSE(RFR)-RFR)

      DEALLOCATE (RM, RMR, F, RF, RFR)
   END FUNCTION fold_derivatives

! **************************************************************************************************
!> \brief Helper routine for diagonalizing anti symmetric matrices
!> \param matrix ...
!> \param evecs ...
!> \param evals ...
! **************************************************************************************************
   SUBROUTINE diag_antisym(matrix, evecs, evals)
      REAL(dp), DIMENSION(:, :)                          :: matrix
      COMPLEX(dp), DIMENSION(:, :)                       :: evecs
      COMPLEX(dp), DIMENSION(:)                          :: evals

      CHARACTER(len=*), PARAMETER :: routineN = 'diag_antisym', routineP = moduleN//':'//routineN

      COMPLEX(dp), DIMENSION(:, :), POINTER              :: matrix_c
      INTEGER                                            :: N
      REAL(dp), DIMENSION(:), POINTER                    :: evals_r

      IF (MAXVAL(ABS(matrix+TRANSPOSE(matrix))) > 1e-14_dp) CPABORT("Expected anti-symmetric matrix")
      N = SIZE(matrix, 1)
      ALLOCATE (matrix_c(N, N), evals_r(N))

      matrix_c = CMPLX(0.0_dp, -matrix, kind=dp)
      CALL zheevd_wrapper(matrix_c, evecs, evals_r)
      evals = CMPLX(0.0_dp, evals_r, kind=dp)

      DEALLOCATE (matrix_c, evals_r)
   END SUBROUTINE diag_antisym

! **************************************************************************************************
!> \brief Helper routine for calling LAPACK zheevd
!> \param matrix ...
!> \param eigenvectors ...
!> \param eigenvalues ...
! **************************************************************************************************
   SUBROUTINE zheevd_wrapper(matrix, eigenvectors, eigenvalues)
      COMPLEX(dp), DIMENSION(:, :)                       :: matrix, eigenvectors
      REAL(dp), DIMENSION(:)                             :: eigenvalues

      CHARACTER(len=*), PARAMETER :: routineN = 'zheevd_wrapper', routineP = moduleN//':'//routineN

      COMPLEX(KIND=dp), DIMENSION(:), POINTER            :: work
      COMPLEX(KIND=dp), DIMENSION(:, :), POINTER         :: A
      INTEGER                                            :: handle, info, liwork, lrwork, lwork, n
      INTEGER, DIMENSION(:), POINTER                     :: iwork
      REAL(KIND=dp), DIMENSION(:), POINTER               :: rwork

      CALL timeset(routineN, handle)

      IF (SIZE(matrix, 1) /= SIZE(matrix, 2)) CPABORT("expected square matrix")
      IF (MAXVAL(ABS(matrix-CONJG(TRANSPOSE(matrix)))) > 1e-14_dp) CPABORT("Expect hermitian matrix")

      n = SIZE(matrix, 1)
      ALLOCATE (iwork(1), rwork(1), work(1), A(n, n))

      A(:, :) = matrix ! ZHEEVD will overwrite A
      ! work space query
      lwork = -1
      lrwork = -1
      liwork = -1

      CALL ZHEEVD('V', 'U', n, A(1, 1), n, eigenvalues(1), &
                  work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)
      lwork = INT(REAL(work(1)))
      lrwork = INT(REAL(rwork(1)))
      liwork = iwork(1)

      DEALLOCATE (iwork, rwork, work)
      ALLOCATE (iwork(liwork))
      iwork(:) = 0
      ALLOCATE (rwork(lrwork))
      rwork(:) = 0.0_dp
      ALLOCATE (work(lwork))
      work(:) = CMPLX(0.0_dp, 0.0_dp, KIND=dp)

      CALL ZHEEVD('V', 'U', n, A(1, 1), n, eigenvalues(1), &
                  work(1), lwork, rwork(1), lrwork, iwork(1), liwork, info)

      eigenvectors = A

      IF (info /= 0) CPABORT("diagonalization failed")

      DEALLOCATE (iwork, rwork, work, A)

      CALL timestop(handle)

   END SUBROUTINE zheevd_wrapper

END MODULE pao_param_direct
