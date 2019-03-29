PROGRAM cp2k

   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CALL main()

CONTAINS

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE main()
      INTEGER, PARAMETER                                 :: m = 5, n = 13

      INTEGER                                            :: a, b, i, j, k
      INTEGER, DIMENSION(:), POINTER                     :: seed
      REAL(dp)                                           :: delta, eps, Gij_num, max_delta
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_A, block_B, block_G, block_Ma, &
                                                            block_Mb, block_N

! number of pao basis functions
! number of primary basis functions

      CALL RANDOM_SEED(size=k)
      ALLOCATE (seed(k))
      seed(:) = 42
      !CALL RANDOM_SEED (PUT=seed)

      ALLOCATE (block_N(n, n), block_A(n, m), block_B(n, m), block_Ma(n, m), block_Mb(n, m), block_G(n, m))
      NULLIFY (block_Ma)

      CALL RANDOM_NUMBER(block_N)
      block_N = 0.5_dp*(block_N+TRANSPOSE(block_N))
      !block_N = 0.0
      !DO i = 1, n
      !   block_N(i,i) = 1.0
      !ENDDO

      CALL RANDOM_NUMBER(block_A)
      !block_A = 0.0
      !DO i = 1, m
      !   block_A(i,i) = 1.0
      !ENDDO
      !block_A(2,1) = 0.2

      max_delta = 0.0_dp
      DO a = 1, n
      DO b = 1, m
         block_Mb = 0.0
         block_Mb(a, b) = 1.0
         CALL calc_AB(block_N, block_A, block_B, block_Mb, block_G)

         eps = 1e-6_dp

         DO i = 1, n
         DO j = 1, m
            !SELECT CASE (pao%num_grad_order)
            !CASE (2) ! calculate derivative to 2th order
            !   Gij_num = -eval_point(block_X, i, j, -eps, pao, ls_scf_env, qs_env)
            !   Gij_num = Gij_num+eval_point(block_X, i, j, +eps, pao, ls_scf_env, qs_env)
            !   Gij_num = Gij_num/(2.0_dp*eps)
            !
            !CASE (4) ! calculate derivative to 4th order
            !   Gij_num = eval_point(block_X, i, j, -2_dp*eps, pao, ls_scf_env, qs_env)
            !   Gij_num = Gij_num-8_dp*eval_point(block_X, i, j, -1_dp*eps, pao, ls_scf_env, qs_env)
            !   Gij_num = Gij_num+8_dp*eval_point(block_X, i, j, +1_dp*eps, pao, ls_scf_env, qs_env)
            !   Gij_num = Gij_num-eval_point(block_X, i, j, +2_dp*eps, pao, ls_scf_env, qs_env)
            !   Gij_num = Gij_num/(12.0_dp*eps)
            !
            !CASE (6) ! calculate derivative to 6th order
            Gij_num = -1_dp*eval_point(block_N, block_A, a, b, i, j, -3_dp*eps)
            Gij_num = Gij_num+9_dp*eval_point(block_N, block_A, a, b, i, j, -2_dp*eps)
            Gij_num = Gij_num-45_dp*eval_point(block_N, block_A, a, b, i, j, -1_dp*eps)
            Gij_num = Gij_num+45_dp*eval_point(block_N, block_A, a, b, i, j, +1_dp*eps)
            Gij_num = Gij_num-9_dp*eval_point(block_N, block_A, a, b, i, j, +2_dp*eps)
            Gij_num = Gij_num+1_dp*eval_point(block_N, block_A, a, b, i, j, +3_dp*eps)
            Gij_num = Gij_num/(60.0_dp*eps)
            !CASE DEFAULT
            !   CPABORT("Unsupported numerical derivative order: "//cp_to_string(pao%num_grad_order))
            !END SELECT

            delta = ABS(Gij_num-block_G(i, j))
            max_delta = MAX(max_delta, delta)
            WRITE (*, *) "gradient check dB(", a, b, ") / dA(", i, j, ") = ", Gij_num, block_G(i, j), delta
         ENDDO
         ENDDO
      ENDDO
      ENDDO
      WRITE (*, *) "max_delta: ", max_delta

   END SUBROUTINE

! **************************************************************************************************
!> \brief ...
!> \param block_N ...
!> \param block_A ...
!> \param a ...
!> \param b ...
!> \param i ...
!> \param j ...
!> \param eps ...
!> \return ...
! **************************************************************************************************
   FUNCTION eval_point(block_N, block_A, a, b, i, j, eps) RESULT(Bab)
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_N, block_A
      INTEGER                                            :: a, b, i, j
      REAL(dp)                                           :: eps, Bab

      INTEGER                                            :: m, n
      REAL(dp)                                           :: old_Aij
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_B, block_G, block_Mb

      n = SIZE(block_A, 1) ! size of primary basis
      m = SIZE(block_A, 2) ! size of pao basis

      NULLIFY (block_Mb, block_G)
      ALLOCATE (block_B(n, m))
      block_B = 123

      old_Aij = block_A(i, j) ! backup old block_A
      block_A(i, j) = block_A(i, j)+eps ! add pertubation

      ! calculate B
      CALL calc_AB(block_N, block_A, block_B, block_Mb, block_G)

      Bab = block_B(a, b)

      ! restore old block_A
      block_A(i, j) = old_Aij

   END FUNCTION eval_point

! **************************************************************************************************
!> \brief ...
!> \param block_N ...
!> \param block_A ...
!> \param block_B ...
!> \param block_Mb ...
!> \param G ...
! **************************************************************************************************
   SUBROUTINE calc_AB(block_N, block_A, block_B, block_Mb, G)

      REAL(dp), DIMENSION(:, :), POINTER                 :: block_N, block_A, block_B, block_Mb, G

      INTEGER                                            :: i, j, k, m, n
      REAL(dp)                                           :: denom, w
      REAL(dp), DIMENSION(:), POINTER                    :: ANNA_evals
      REAL(dp), DIMENSION(:, :), POINTER                 :: ANNA, ANNA_evecs, ANNA_inv, D, M1, M2, &
                                                            M3, M4, M5, NN

      n = SIZE(block_A, 1) ! size of primary basis
      m = SIZE(block_A, 2) ! size of pao basis

      !! restrain pao basis vectors to unit norm
      !IF (PRESENT(penalty)) THEN
      !   DO i = 1, m
      !      w = 1.0_dp-SUM(block_A(:, i)**2)
      !      penalty = penalty+pao%penalty_strength*w**2
      !   ENDDO
      !ENDIF

      ALLOCATE (NN(n, n), ANNA(m, m))
      NN = MATMUL(block_N, block_N) ! it's actually S^{-1}
      ANNA = MATMUL(MATMUL(TRANSPOSE(block_A), NN), block_A)

      ! diagonalize ANNA
      ALLOCATE (ANNA_evecs(m, m), ANNA_evals(m))
      ANNA_evecs(:, :) = ANNA
      CALL diamat_all(ANNA_evecs, ANNA_evals)
      IF (MINVAL(ABS(ANNA_evals)) < 1e-10_dp) CPABORT("PAO basis singualar.")

      ! build ANNA_inv
      ALLOCATE (ANNA_inv(m, m))
      ANNA_inv(:, :) = 0.0_dp
      DO k = 1, m
         w = 1.0_dp/ANNA_evals(k)
         DO i = 1, m
         DO j = 1, m
            ANNA_inv(i, j) = ANNA_inv(i, j)+w*ANNA_evecs(i, k)*ANNA_evecs(j, k)
         ENDDO
         ENDDO
      ENDDO

      !B = 1/S * A * 1/(A^T 1/S A)
      block_B = MATMUL(MATMUL(NN, block_A), ANNA_inv)

      ! TURNING POINT (if calc grad) ------------------------------------------
      IF (ASSOCIATED(G)) THEN
         WRITE (*, *) "ANNA_evals: ", ANNA_evals

         G(:, :) = 0.0_dp

         !IF (PRESENT(penalty)) THEN
         !   DO i = 1, m
         !      w = 1.0_dp-SUM(block_A(:, i)**2)
         !      G(:, i) = -4.0_dp*pao%penalty_strength*w*block_A(:, i)
         !   ENDDO
         !ENDIF

         !IF (ASSOCIATED(block_Ma)) THEN
         !   G = G+block_Ma
         !ENDIF

         IF (ASSOCIATED(block_Mb)) THEN
            G = G+MATMUL(MATMUL(NN, block_Mb), ANNA_inv)

            ! calculate derivatives dAA_inv/ dAA
            ALLOCATE (D(m, m), M1(m, m), M2(m, m), M3(m, m), M4(m, m), M5(m, m))

            DO i = 1, m
            DO j = 1, m
               denom = ANNA_evals(i)-ANNA_evals(j)
               IF (i == j) THEN
                  D(i, i) = -1.0_dp/ANNA_evals(i)**2 ! diagonal elements
               ELSE IF (ABS(denom) > 1e-10_dp) THEN
                  D(i, j) = (1.0_dp/ANNA_evals(i)-1.0_dp/ANNA_evals(j))/denom
               ELSE
                  D(i, j) = -1.0_dp ! limit according to L'Hospital's rule
               ENDIF
            ENDDO
            ENDDO

            M1 = MATMUL(MATMUL(TRANSPOSE(block_A), NN), block_Mb)
            M2 = MATMUL(MATMUL(TRANSPOSE(ANNA_evecs), M1), ANNA_evecs)
            M3 = M2*D ! Hadamard product
            M4 = MATMUL(MATMUL(ANNA_evecs, M3), TRANSPOSE(ANNA_evecs))
            M5 = 0.5_dp*(M4+TRANSPOSE(M4))
            G = G+2.0_dp*MATMUL(MATMUL(NN, block_A), M5)

            DEALLOCATE (D, M1, M2, M3, M4, M5)
         ENDIF
      ENDIF
   END SUBROUTINE calc_AB

END PROGRAM
