!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Gradient of LNV (Li, Nunes, Vanderbilt) energy
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_lnv
   USE cp_control_types,                ONLY: dft_control_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_create, dbcsr_get_info, dbcsr_multiply, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_diag_blocks, dbcsr_scale, dbcsr_type
   USE dm_ls_scf_qs,                    ONLY: matrix_decluster
   USE dm_ls_scf_types,                 ONLY: ls_mstruct_type,&
                                              ls_scf_env_type
   USE kinds,                           ONLY: dp
   USE pao_types,                       ONLY: pao_env_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_lnv'

   PUBLIC :: pao_calc_grad_lnv_wrt_AB

CONTAINS

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

END MODULE pao_lnv
