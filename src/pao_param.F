!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Front-End for any PAO parametrization
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_to_string
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, dbcsr_frobenius_norm, &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_multiply, &
        dbcsr_p_type, dbcsr_release, dbcsr_reserve_diag_blocks, dbcsr_scale, dbcsr_type
   USE dm_ls_scf_qs,                    ONLY: matrix_decluster
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_max
   USE pao_input,                       ONLY: pao_exp_param,&
                                              pao_fock_param,&
                                              pao_gth_param,&
                                              pao_rotinv_param
   USE pao_param_exp,                   ONLY: pao_calc_U_exp,&
                                              pao_param_count_exp,&
                                              pao_param_finalize_exp,&
                                              pao_param_init_exp,&
                                              pao_param_initguess_exp
   USE pao_param_gth,                   ONLY: pao_calc_U_gth,&
                                              pao_param_count_gth,&
                                              pao_param_finalize_gth,&
                                              pao_param_init_gth,&
                                              pao_param_initguess_gth
   USE pao_param_linpot,                ONLY: pao_calc_U_linpot,&
                                              pao_param_count_linpot,&
                                              pao_param_finalize_linpot,&
                                              pao_param_init_linpot,&
                                              pao_param_initguess_linpot
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_param'

   PUBLIC :: pao_update_AB, pao_param_count, pao_param_initial_guess
   PUBLIC :: pao_param_init, pao_param_finalize, pao_calc_U

CONTAINS

! **************************************************************************************************
!> \brief Takes current matrix_X and recalculates derived matrices U, A, and B.
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param penalty ...
! **************************************************************************************************
   SUBROUTINE pao_update_AB(pao, qs_env, ls_scf_env, penalty)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      REAL(dp), INTENT(OUT), OPTIONAL                    :: penalty

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_update_AB'

      INTEGER                                            :: acol, arow, handle, iatom
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_N, block_N_inv, &
                                                            block_U, block_Y
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: matrix_U
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct

      CALL pao_calc_U(pao, qs_env, ls_scf_env, gradient=.FALSE., penalty=penalty) !update matrix_U = Function of matrix_X

      ! pao%matrix_U uses pao%diag_distribution, need to redistribute using distribution of matrix_s
      CALL get_qs_env(qs_env, matrix_s=matrix_s)
      CALL dbcsr_create(matrix_U, matrix_type="N", template=matrix_s(1)%matrix)
      CALL dbcsr_reserve_diag_blocks(matrix_U)
      CALL dbcsr_complete_redistribute(pao%matrix_U, matrix_U)

      ! Multiplying diagonal matrices is a local operation.
      ! To take advantage of this we're using an iterator instead of calling dbcsr_multiply().
!$OMP PARALLEL DEFAULT(NONE) SHARED(pao,ls_mstruct,matrix_U) &
!$OMP PRIVATE(iter,arow,acol,iatom,block_U,block_Y,block_A,block_B,block_N,block_N_inv,found)
      CALL dbcsr_iterator_start(iter, matrix_U)
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

      CALL dbcsr_release(matrix_U)

      CALL timestop(handle)
   END SUBROUTINE pao_update_AB

! **************************************************************************************************
!> \brief Initialize PAO parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_init_exp(pao, qs_env)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_init_linpot(pao, qs_env)
      CASE (pao_gth_param)
         CALL pao_param_init_gth(pao, qs_env)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_init

! **************************************************************************************************
!> \brief Finalize PAO parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_finalize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_finalize_exp(pao)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_finalize_linpot(pao)
      CASE (pao_gth_param)
         CALL pao_param_finalize_gth(pao)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_finalize

! **************************************************************************************************
!> \brief Returns the number of parameters for given atomic kind
!> \param pao ...
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count(pao, qs_env, ikind, nparams)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_count', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_count_exp(qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_count_linpot(pao, qs_env, ikind=ikind, nparams=nparams)
      CASE (pao_gth_param)
         CALL pao_param_count_gth(qs_env, ikind=ikind, nparams=nparams)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pao_param_count

! **************************************************************************************************
!> \brief Fills matrix_X with an initial guess
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_initial_guess(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initial_guess', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      REAL(dp)                                           :: norm
      TYPE(dbcsr_type)                                   :: matrix_tmp

      CALL timeset(routineN, handle)

      SELECT CASE (pao%parameterization)
      CASE (pao_exp_param)
         CALL pao_param_initguess_exp(pao)
      CASE (pao_fock_param, pao_rotinv_param)
         CALL pao_param_initguess_linpot(pao, qs_env)
      CASE (pao_gth_param)
         CALL pao_param_initguess_gth(pao)
      CASE DEFAULT
         CPABORT("PAO: unkown parametrization")
      END SELECT

      norm = dbcsr_frobenius_norm(pao%matrix_X)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Made initial guess for matrix_X with norm:", norm

      IF (pao%precondition) THEN
         !TODO: multiplying a matrix into itself while retaining sparsity seems to be broken
         CALL dbcsr_copy(matrix_tmp, pao%matrix_X)
         CALL dbcsr_multiply("N", "N", 1.0_dp, pao%matrix_precon, matrix_tmp, &
                             0.0_dp, pao%matrix_X, retain_sparsity=.TRUE.)
         CALL dbcsr_release(matrix_tmp)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE pao_param_initial_guess

! **************************************************************************************************
!> \brief Calculate new matrix U
!> \param pao ...
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_G ...
!> \param penalty ...
!> \param forces ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U(pao, qs_env, ls_scf_env, gradient, penalty, forces)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type)                              :: ls_scf_env
      LOGICAL                                            :: gradient
      REAL(dp), INTENT(OUT), OPTIONAL                    :: penalty
      REAL(dp), DIMENSION(:, :), INTENT(INOUT), OPTIONAL :: forces

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U'

      INTEGER                                            :: handle
      TYPE(dbcsr_type)                                   :: matrix_M

      CALL timeset(routineN, handle)

      IF (PRESENT(penalty)) penalty = 0.0_dp

      IF (gradient) THEN
         CALL pao_calc_grad_lnv_wrt_U(qs_env, ls_scf_env, matrix_M)

         SELECT CASE (pao%parameterization)
         CASE (pao_exp_param)
            CALL pao_calc_U_exp(pao, matrix_M, pao%matrix_G)
         CASE (pao_fock_param, pao_rotinv_param)
            CALL pao_calc_U_linpot(pao, qs_env, penalty, matrix_M, pao%matrix_G, forces)
         CASE (pao_gth_param)
            CALL pao_calc_U_gth(pao, penalty, matrix_M, pao%matrix_G)
         CASE DEFAULT
            CPABORT("PAO: unkown parametrization")
         END SELECT

         CALL dbcsr_release(matrix_M)

      ELSE

         SELECT CASE (pao%parameterization)
         CASE (pao_exp_param)
            CALL pao_calc_U_exp(pao)
         CASE (pao_fock_param, pao_rotinv_param)
            CALL pao_calc_U_linpot(pao, qs_env, penalty)
         CASE (pao_gth_param)
            CALL pao_calc_U_gth(pao, penalty)
         CASE DEFAULT
            CPABORT("PAO: unkown parametrization")
         END SELECT

      ENDIF

      CALL pao_assert_unitary(pao)
      CALL timestop(handle)
   END SUBROUTINE pao_calc_U

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
!> \brief Helper routine, calculates partial derivative dE/dA and dE/dB
!> \param qs_env ...
!> \param ls_scf_env ...
!> \param matrix_Ma the derivate wrt A, matrix uses s_matrix-distribution.
!> \param matrix_Mb the derivate wrt B, matrix uses s_matrix-distribution.
! **************************************************************************************************
   SUBROUTINE pao_calc_grad_lnv_wrt_AB(qs_env, ls_scf_env, matrix_Ma, matrix_Mb)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(ls_scf_env_type), TARGET                      :: ls_scf_env
      TYPE(dbcsr_type)                                   :: matrix_Ma, matrix_Mb

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_calc_grad_lnv_wrt_AB', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nspin
      INTEGER, DIMENSION(:), POINTER                     :: pao_blk_sizes
      REAL(KIND=dp)                                      :: filter_eps
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, rho_ao
      TYPE(dbcsr_type) :: matrix_HB, matrix_HPS, matrix_M, matrix_M1, matrix_M1_dc, matrix_M2, &
         matrix_M2_dc, matrix_M3, matrix_M3_dc, matrix_PA, matrix_PH, matrix_PHP, matrix_PSP, &
         matrix_SB, matrix_SP
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(ls_mstruct_type), POINTER                     :: ls_mstruct
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_rho_type), POINTER                         :: rho

      CALL timeset(routineN, handle)

      ls_mstruct => ls_scf_env%ls_mstruct
      pao => ls_scf_env%pao_env

      CALL get_qs_env(qs_env, &
                      rho=rho, &
                      matrix_ks=matrix_ks, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control)
      CALL qs_rho_get(rho, rho_ao=rho_ao)
      nspin = dft_control%nspins
      filter_eps = ls_scf_env%eps_filter

      CALL dbcsr_get_info(ls_mstruct%matrix_A, col_blk_size=pao_blk_sizes)

      IF (nspin /= 1) CPABORT("open shell not yet implemented")
      !TODO: handle openshell case properly

      ! Notation according to equation (4.6) on page 50 from:
      ! https://dx.doi.org/10.3929%2Fethz-a-010819495

      !---------------------------------------------------------------------------
      ! calculate need products in pao basis
      CALL dbcsr_create(matrix_PH, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_ks(1), &
                          0.0_dp, matrix_PH, filter_eps=filter_eps)

      CALL dbcsr_create(matrix_PHP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_PH, ls_scf_env%matrix_p(1), &
                          0.0_dp, matrix_PHP, filter_eps=filter_eps)

      CALL dbcsr_create(matrix_SP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s, ls_scf_env%matrix_p(1), &
                          0.0_dp, matrix_SP, filter_eps=filter_eps)

      IF (nspin == 1) CALL dbcsr_scale(matrix_SP, 0.5_dp)

      CALL dbcsr_create(matrix_HPS, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "T", 1.0_dp, ls_scf_env%matrix_ks(1), matrix_SP, &
                          0.0_dp, matrix_HPS, filter_eps=filter_eps)

      CALL dbcsr_create(matrix_PSP, template=ls_scf_env%matrix_s, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), matrix_SP, &
                          0.0_dp, matrix_PSP, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! M1 = dE_lnv / dP_pao
      CALL dbcsr_create(matrix_M1, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL dbcsr_multiply("N", "T", 3.0_dp, ls_scf_env%matrix_ks(1), matrix_SP, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", 3.0_dp, matrix_SP, ls_scf_env%matrix_ks(1), &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", -2.0_dp, matrix_HPS, matrix_SP, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_SP, matrix_HPS, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", -2.0_dp, matrix_SP, matrix_HPS, &
                          1.0_dp, matrix_M1, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL dbcsr_create(matrix_M1_dc, &
                        template=matrix_s(1)%matrix, &
                        row_blk_size=pao_blk_sizes, &
                        col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M1_dc, matrix_M1, ls_mstruct)

      !---------------------------------------------------------------------------
      ! M2 = dE_lnv / dH
      CALL dbcsr_create(matrix_M2, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL dbcsr_add(matrix_M2, matrix_PSP, 1.0_dp, 3.0_dp)

      CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_PSP, matrix_SP, &
                          1.0_dp, matrix_M2, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL dbcsr_create(matrix_M2_dc, &
                        template=matrix_s(1)%matrix, &
                        row_blk_size=pao_blk_sizes, &
                        col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M2_dc, matrix_M2, ls_mstruct)

      !---------------------------------------------------------------------------
      ! M3 = dE_lnv / dS
      CALL dbcsr_create(matrix_M3, template=ls_scf_env%matrix_s, matrix_type="N")

      CALL dbcsr_add(matrix_M3, matrix_PHP, 1.0_dp, 3.0_dp)

      CALL dbcsr_multiply("N", "N", -2.0_dp, matrix_PHP, matrix_SP, &
                          1.0_dp, matrix_M3, filter_eps=filter_eps)

      CALL dbcsr_multiply("N", "T", -2.0_dp, matrix_PSP, matrix_PH, &
                          1.0_dp, matrix_M3, filter_eps=filter_eps)

      ! reverse possible molecular clustering
      CALL dbcsr_create(matrix_M3_dc, &
                        template=matrix_s(1)%matrix, &
                        row_blk_size=pao_blk_sizes, &
                        col_blk_size=pao_blk_sizes)
      CALL matrix_decluster(matrix_M3_dc, matrix_M3, ls_mstruct)

      !---------------------------------------------------------------------------
      ! assemble Ma and Mb
      ! matrix_Ma = dE_lnv / dA = P * A * M1
      ! matrix_Mb = dE_lnv / dB = H * B * M2  +  S * B * M3
      CALL dbcsr_create(matrix_Ma, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL dbcsr_reserve_diag_blocks(matrix_Ma)
      CALL dbcsr_create(matrix_Mb, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL dbcsr_reserve_diag_blocks(matrix_Mb)

      !---------------------------------------------------------------------------
      ! combine M1 with matrices from primary basis
      CALL dbcsr_create(matrix_PA, template=ls_mstruct%matrix_A, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, rho_ao(1)%matrix, ls_mstruct%matrix_A, &
                          0.0_dp, matrix_PA, filter_eps=filter_eps)

      ! matrix_Ma = P * A * M1
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_PA, matrix_M1_dc, &
                          0.0_dp, matrix_Ma, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! combine M2 with matrices from primary basis
      CALL dbcsr_create(matrix_HB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_ks(1)%matrix, ls_mstruct%matrix_B, &
                          0.0_dp, matrix_HB, filter_eps=filter_eps)

      ! matrix_Mb = H * B * M2
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_HB, matrix_M2_dc, &
                          0.0_dp, matrix_Mb, filter_eps=filter_eps)

      !---------------------------------------------------------------------------
      ! combine M3 with matrices from primary basis
      CALL dbcsr_create(matrix_SB, template=ls_mstruct%matrix_B, matrix_type="N")
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_s(1)%matrix, ls_mstruct%matrix_B, &
                          0.0_dp, matrix_SB, filter_eps=filter_eps)

      IF (nspin == 1) CALL dbcsr_scale(matrix_SB, 0.5_dp)

      ! matrix_Mb += S * B * M3
      CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_SB, matrix_M3_dc, &
                          1.0_dp, matrix_Mb, filter_eps=filter_eps)

      IF (nspin == 1) CALL dbcsr_scale(matrix_Ma, 2.0_dp)
      IF (nspin == 1) CALL dbcsr_scale(matrix_Mb, 2.0_dp)

      !---------------------------------------------------------------------------
      ! cleanup: TODO release matrices as early as possible
      CALL dbcsr_release(matrix_PH)
      CALL dbcsr_release(matrix_PHP)
      CALL dbcsr_release(matrix_SP)
      CALL dbcsr_release(matrix_HPS)
      CALL dbcsr_release(matrix_PSP)
      CALL dbcsr_release(matrix_M)
      CALL dbcsr_release(matrix_M1)
      CALL dbcsr_release(matrix_M2)
      CALL dbcsr_release(matrix_M3)
      CALL dbcsr_release(matrix_M1_dc)
      CALL dbcsr_release(matrix_M2_dc)
      CALL dbcsr_release(matrix_M3_dc)
      CALL dbcsr_release(matrix_PA)
      CALL dbcsr_release(matrix_HB)
      CALL dbcsr_release(matrix_SB)

      CALL timestop(handle)
   END SUBROUTINE pao_calc_grad_lnv_wrt_AB

! **************************************************************************************************
!> \brief Debugging routine, check unitaryness of U
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_assert_unitary(pao)
      TYPE(pao_env_type), POINTER                        :: pao

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
      CALL dbcsr_iterator_start(iter, pao%matrix_U)
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

      CALL dbcsr_get_info(pao%matrix_U, group=group)
      CALL mp_max(delta_max, group)
      IF (pao%iw > 0) WRITE (pao%iw, *) 'PAO| checked unitaryness, max delta:', delta_max
      IF (delta_max > pao%check_unitary_tol) &
         CPABORT("Found bad unitaryness:"//cp_to_string(delta_max))

      CALL timestop(handle)
   END SUBROUTINE pao_assert_unitary

END MODULE pao_param
