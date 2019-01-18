!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common routines for PAO parametrizations that rely on the unitary matrix U.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_unitary
   USE cp_log_handling,                 ONLY: cp_to_string
   USE dbcsr_api,                       ONLY: &
        dbcsr_complete_redistribute, dbcsr_create, dbcsr_get_block_p, dbcsr_get_info, &
        dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, dbcsr_iterator_start, &
        dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_diag_blocks, dbcsr_type
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_max
   USE pao_input,                       ONLY: pao_exp_param,&
                                              pao_fock_param,&
                                              pao_gth_param,&
                                              pao_rotinv_param
   USE pao_lnv,                         ONLY: pao_calc_grad_lnv_wrt_AB
   USE pao_param_exp,                   ONLY: pao_calc_U_exp
   USE pao_param_gth,                   ONLY: pao_calc_U_gth
   USE pao_param_linpot,                ONLY: pao_calc_U_linpot
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param_unitary'

   PUBLIC :: pao_calc_AB_unitary

CONTAINS

! **************************************************************************************************
!> \brief Takes current matrix_X and calculates the matrices A and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param gradient ...
!> \param penalty ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_calc_AB_unitary(pao, qs_env, ls_scf_env, gradient, penalty, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      LOGICAL, INTENT(IN)                                :: gradient
      REAL(dp), INTENT(INOUT), OPTIONAL                  :: penalty
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_AB_unitary'

      INTEGER                                            :: acol, arow, handle, iatom
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_N, block_N_inv, &
                                                            block_U, block_Y
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_M, matrix_U, matrix_U_sdist
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct

      CALL timeset(routineN, handle)
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      ls_mstruct => ls_scf_env%ls_mstruct

      ! --------------------------------------------------------------------------------------------
      ! calculate matrix_U using diag_distribution
      CALL dbcsr_create(matrix_U, matrix_type="N", dist=pao%diag_distribution, template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U)

      ! TODO: clean up this duplicationed SELECT block.
      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_U(qs_env, ls_scf_env, matrix_M)

         SELECT CASE (pao%parameterization)
         CASE (pao_exp_param)
            CALL pao_calc_U_exp(pao, matrix_U, matrix_M, pao%matrix_G)
         CASE (pao_fock_param, pao_rotinv_param)
            CALL pao_calc_U_linpot(pao, qs_env, matrix_U, matrix_M, pao%matrix_G, penalty, forces)
         CASE (pao_gth_param)
            CALL pao_calc_U_gth(pao, matrix_U, matrix_M, pao%matrix_G, penalty)
         CASE DEFAULT
            CPABORT("PAO: unkown parametrization")
         END SELECT

         CALL dbcsr_release(matrix_M)

      ELSE

         SELECT CASE (pao%parameterization)
         CASE (pao_exp_param)
            CALL pao_calc_U_exp(pao, matrix_U)
         CASE (pao_fock_param, pao_rotinv_param)
            CALL pao_calc_U_linpot(pao, qs_env, matrix_U, penalty=penalty)
         CASE (pao_gth_param)
            CALL pao_calc_U_gth(pao, matrix_U, penalty=penalty)
         CASE DEFAULT
            CPABORT("PAO: unkown parametrization")
         END SELECT

      ENDIF

      !TODO: move the remainder of this into its own subroutine
      !TODO: replace this routine with pao_calc_AB_(exp/linpot/gth)

      CALL pao_assert_unitary(pao, matrix_U) ! sanity check matrix U

      ! --------------------------------------------------------------------------------------------
      ! redistribute matrix_U from diag_distribution to distribution of matrix_s
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_create(matrix_U_sdist, matrix_type="N", template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U_sdist)
      CALL dbcsr_complete_redistribute(matrix_U, matrix_U_sdist)
      CALL dbcsr_release(matrix_U)

      ! --------------------------------------------------------------------------------------------
      ! calculate matrix A and B from matrix U
      ! Multiplying diagonal matrices is a local operation.
      ! To take advantage of this we're using an iterator instead of calling dbcsr_multiply().
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,matrix_U_sdist) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_U,block_Y,block_A,block_B,block_N,block_N_inv,found)
      CALL dbcsr_iterator_start(iter, matrix_U_sdist)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_U)
         iatom = arow; CPASSERT(arow == acol)

         CALL dbcsr_get_block_p(matrix=pao%matrix_Y, row=iatom, col=iatom, block=block_Y, found=found)
         CPASSERT(ASSOCIATED(block_Y))

         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_A, row=iatom, col=iatom, block=block_A, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N_inv, row=iatom, col=iatom, block=block_N_inv, found=found)
         CPASSERT(ASSOCIATED(block_A) .AND. ASSOCIATED(block_N_inv))

         CALL dbcsr_get_block_p(matrix=ls_mstruct%matrix_B, row=iatom, col=iatom, block=block_B, found=found)
         CALL dbcsr_get_block_p(matrix=pao%matrix_N, row=iatom, col=iatom, block=block_N, found=found)
         CPASSERT(ASSOCIATED(block_B) .AND. ASSOCIATED(block_N))

         block_A = MATMUL(MATMUL(block_N_inv, block_U), block_Y)
         block_B = MATMUL(MATMUL(block_N, block_U), block_Y)
      ENDDO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_release(matrix_U_sdist)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_AB_unitary

! **************************************************************************************************
!> \brief Helper routine, calculates partial derivative dE/dU
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_M_diag the derivate wrt U, matrix uses pao%diag_distribution
! **************************************************************************************************
   SUBROUTINE pao_calc_grad_lnv_wrt_U(qs_env, ls_scf_env, matrix_M_diag)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(dbcsr_type)                                   :: matrix_M_diag

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad_lnv_wrt_U', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(KIND=dp)                                      :: filter_eps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_M, matrix_Ma, matrix_Mb, matrix_NM
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env
      filter_eps = ls_scf_env%eps_filter
      CALL get_qs_env(qs_env, matrix_s=matrix_s)

      CALL pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)

      ! Calculation uses distr. of matrix_s, afterwards we redistribute to pao%diag_distribution.
      CALL dbcsr_create(matrix_M, template=matrix_s(1)%matrix, matrix_type="N")
      CALL dbcsr_reserve_diag_blocks(matrix_M)

      CALL dbcsr_create(matrix_NM, template=ls_mstruct%matrix_A, matrix_type="N")

      CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N_inv, matrix_Ma, &
                          1.0_dp, matrix_NM, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_N, matrix_Mb, &
                          1.0_dp, matrix_NM, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", 1.0_dp, matrix_NM, pao%matrix_Y, &
                          1.0_dp, matrix_M, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! redistribute using pao%diag_distribution
      CALL dbcsr_create(matrix_M_diag, &
                        name="PAO matrix_M", &
                        matrix_type="N", &
                        dist=pao%diag_distribution, &
                        template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_M_diag)
      CALL dbcsr_complete_redistribute(matrix_M, matrix_M_diag)

      !---------------------------------------------------------------------------
      ! cleanup:
      CALL dbcsr_release(matrix_M)
      CALL dbcsr_release(matrix_Ma)
      CALL dbcsr_release(matrix_Mb)
      CALL dbcsr_release(matrix_NM)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_grad_lnv_wrt_U

! **************************************************************************************************
!> \brief Debugging routine, check unitaryness of U
!> \param pao ...
!> \param matrix_U ...
! **************************************************************************************************
   SUBROUTINE pao_assert_unitary(pao, matrix_U)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(dbcsr_type)                                   :: matrix_U

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_assert_unitary', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: acol, arow, group, handle, i, iatom, M, N
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pao, blk_sizes_pri
      REAL(dp)                                           :: delta_max
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_test, tmp1, tmp2
      TYPE(dbcsr_iterator_type)                          :: iter

      IF (pao%check_unitary_tol < 0.0_dp) RETURN ! no checking

      CALL timeset(routineN, handle)
      delta_max = 0.0_dp

      CALL dbcsr_get_info(pao%matrix_Y, row_blk_size=blk_sizes_pri, col_blk_size=blk_sizes_pao)

!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,blk_sizes_pri,blk_sizes_pao,delta_max) &
!$OMP PRIVATE(iter,arow,acol,iatom,N,M,block_test,tmp1,tmp2)
      CALL dbcsr_iterator_start(iter, matrix_U)
      DO WHILE (dbcsr_iterator_blocks_left(iter))
         CALL dbcsr_iterator_next_block(iter, arow, acol, block_test)
         iatom = arow; CPASSERT(arow == acol)
         N = blk_sizes_pri(iatom) ! size of primary basis
         M = blk_sizes_pao(iatom) ! size of pao basis

         ! we only need the upper left "PAO-corner" to be unitary
         ALLOCATE (tmp1(N, M), tmp2(M, M))
         tmp1 = block_test(:, 1:M)
         tmp2 = MATMUL(TRANSPOSE(tmp1), tmp1)
         DO i = 1, M
            tmp2(i, i) = tmp2(i, i)-1.0_dp
         ENDDO

!$OMP ATOMIC
         delta_max = MAX(delta_max, MAXVAL(ABS(tmp2)))

         DEALLOCATE (tmp1, tmp2)
      END DO
      CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

      CALL dbcsr_get_info(matrix_U, group=group)
      CALL mp_max(delta_max, group)
      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO| checked unitaryness, max delta:', delta_max
      IF (delta_max > pao%check_unitary_tol) &
         CPABORT("Found bad unitaryness:"//cp_to_string(delta_max))

      CALL timestop(handle)
   END SUBROUTINE pao_assert_unitary

END MODULE pao_param_unitary
