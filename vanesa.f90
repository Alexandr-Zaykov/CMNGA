! ==============================================================================
! MNACMEA: Multi-Niche Adaptive Covariance Matrix Evolution Algorithm
! ==============================================================================
! FIXED VERSION - Addresses memory corruption issues
! ==============================================================================

MODULE blas_external
    IMPLICIT NONE
    EXTERNAL :: sgemv, sger, sgetrf, sgetri, sgecon, spotrf, spotri
END MODULE blas_external

MODULE fitness_functions
    IMPLICIT NONE
    
CONTAINS
    
    FUNCTION fitness_five_peaks(x, dims) RESULT(f)
        REAL, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        REAL :: f
        REAL :: a(5), b(5), h(5), w(5)
        INTEGER :: j
        
        DATA a /-20.0, -5.0, 0.0, 30.0, 30.0/
        DATA b /-20.0, -25.0, 30.0, 0.0, -30.0/
        DATA h /0.4, 0.2, 0.7, 1.0, 0.05/
        DATA w /0.02, 0.5, 0.01, 2.0, 0.1/
        
        IF (dims == 1) THEN
            f = SIN(5.0 * x(1))**6 * EXP(-x(1)**2) + 1.0
        ELSE
            f = 0.0
            DO j = 1, 5
                f = f + h(j) / (1.0 + w(j) * ((x(1)-a(j))**2 + (x(2)-b(j))**2))
            END DO
            f = f + 1.0
            
            IF (dims > 2) THEN
                DO j = 3, dims
                    f = f * EXP(-(x(j)**2)/4.0)
                END DO
            END IF
        END IF
    END FUNCTION fitness_five_peaks
    
    FUNCTION fitness_rastrigin(x, dims) RESULT(f)
        REAL, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        REAL :: f
        INTEGER :: i
        REAL, PARAMETER :: pi = 3.14159265359
        
        f = 10.0 * dims
        DO i = 1, dims
            f = f + x(i)**2 - 10.0 * COS(2.0 * pi * x(i))
        END DO
        f = -f + 100.0
    END FUNCTION fitness_rastrigin
    
    FUNCTION fitness_sphere(x, dims) RESULT(f)
        REAL, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        REAL :: f
        
        f = -SUM(x**2) + 100.0
    END FUNCTION fitness_sphere
    
    FUNCTION fitness_schwefel(x, dims) RESULT(f)
        REAL, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        REAL :: f
        INTEGER :: i
        
        f = 0.0
        DO i = 1, dims
            f = f + x(i) * SIN(SQRT(ABS(x(i))))
        END DO
        f = 418.9829 * dims + f
    END FUNCTION fitness_schwefel
    
    FUNCTION fitness_equal_peaks(x, dims) RESULT(f)
        REAL, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        REAL :: f
        REAL :: peaks(4, 2)
        INTEGER :: i
        REAL :: dist, peak_val
        
        DATA peaks /-20.0, -20.0, -20.0, 20.0, 20.0, -20.0, 20.0, 20.0/
        
        f = 0.0
        DO i = 1, 4
            dist = SQRT((x(1) - peaks(i,1))**2 + (x(2) - peaks(i,2))**2)
            peak_val = EXP(-dist**2 / 50.0)
            f = MAX(f, peak_val)
        END DO
        f = f + 1.0
    END FUNCTION fitness_equal_peaks
    
END MODULE fitness_functions

! ==============================================================================

MODULE bfgs_refinement
    USE blas_external
    IMPLICIT NONE
    
CONTAINS

    SUBROUTINE refine_niche_bfgs(x_start, b_inv_init, fitness_val, num_vars, &
                                 fitness_func_id, lower, upper, trust_radius, num_eval)
        REAL, DIMENSION(num_vars), INTENT(INOUT) :: x_start
        REAL, DIMENSION(num_vars, num_vars), INTENT(IN) :: b_inv_init
        REAL, INTENT(INOUT) :: fitness_val
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        REAL, INTENT(IN) :: lower, upper, trust_radius
        INTEGER, INTENT(INOUT) :: num_eval
        
        INTEGER, PARAMETER :: max_iter = 50
        REAL, PARAMETER :: grad_tol = 1e-5
        REAL, PARAMETER :: step_tol = 1e-8
        REAL, PARAMETER :: f_tol = 1e-8
        REAL, PARAMETER :: min_grad_for_refinement = 1e-4
        
        REAL, DIMENSION(num_vars) :: x, x_new, grad, grad_new, search_dir, s, y, x_origin
        ! FIX 1: Use proper contiguous arrays for BLAS temporaries
        REAL, DIMENSION(num_vars) :: temp_vec1
        REAL, DIMENSION(num_vars, num_vars) :: b_inv, b_inv_new
        REAL :: f, f_new, alpha, grad_norm, s_dot_y, rho, y_dot_b_inv_y
        REAL :: improvement, dist_from_origin, f_origin
        INTEGER :: iter, i, j
        LOGICAL :: converged
        
        ! Explicit copy to avoid array slice issues
        DO i = 1, num_vars
            x(i) = x_start(i)
            x_origin(i) = x_start(i)
        END DO
        f = fitness_val
        f_origin = fitness_val
        
        ! FIX 2: Explicit copy of matrix to ensure contiguous memory
        DO i = 1, num_vars
            DO j = 1, num_vars
                b_inv(i,j) = b_inv_init(i,j)
            END DO
        END DO
        
        converged = .FALSE.
        
        CALL compute_gradient_fd(x, grad, num_vars, fitness_func_id, lower, upper, num_eval)
        DO i = 1, num_vars
            grad(i) = -grad(i)
        END DO
        
        grad_norm = SQRT(DOT_PRODUCT(grad, grad))
        
        PRINT *, "  Starting BFGS refinement at x =", x
        PRINT *, "  Initial fitness =", f
        PRINT *, "  Initial gradient norm =", grad_norm
        PRINT *, "  Trust region radius =", trust_radius
        
        IF (grad_norm < min_grad_for_refinement) THEN
            PRINT *, "  Skipping: Already converged (gradient too small)"
            PRINT *
            RETURN
        END IF
        
        DO iter = 1, max_iter
            grad_norm = SQRT(DOT_PRODUCT(grad, grad))
            
            IF (grad_norm < grad_tol) THEN
                PRINT *, "  BFGS converged: gradient tolerance reached"
                converged = .TRUE.
                EXIT
            END IF
            
            ! FIX 3: Use contiguous temp_vec1 for search direction
            CALL sgemv('N', num_vars, num_vars, -1.0, b_inv, num_vars, grad, 1, 0.0, search_dir, 1)
            
            CALL line_search_backtracking(x, search_dir, f, grad, alpha, x_new, f_new, &
                                         num_vars, fitness_func_id, lower, upper, num_eval)
            
            dist_from_origin = 0.0
            DO i = 1, num_vars
                dist_from_origin = dist_from_origin + (x_new(i) - x_origin(i))**2
            END DO
            dist_from_origin = SQRT(dist_from_origin)
            
            IF (dist_from_origin > trust_radius) THEN
                PRINT *, "  BFGS stopped: Exceeded trust region (dist =", dist_from_origin, ")"
                PRINT *, "  This prevents jumping to different basin - keeping previous best"
                converged = .TRUE.
                EXIT
            END IF
            
            improvement = f_new - f
            
            IF (improvement < -1e-6 .AND. iter > 3) THEN
                PRINT *, "  BFGS stopped: Fitness degrading (change =", improvement, ")"
                converged = .TRUE.
                EXIT
            END IF
            
            IF (ABS(improvement) < f_tol .AND. iter > 5) THEN
                PRINT *, "  BFGS converged: fitness change too small"
                converged = .TRUE.
                DO i = 1, num_vars
                    x(i) = x_new(i)
                END DO
                f = f_new
                EXIT
            END IF
            
            CALL compute_gradient_fd(x_new, grad_new, num_vars, fitness_func_id, lower, upper, num_eval)
            DO i = 1, num_vars
                grad_new(i) = -grad_new(i)
            END DO
            
            DO i = 1, num_vars
                s(i) = x_new(i) - x(i)
                y(i) = grad_new(i) - grad(i)
            END DO
            
            IF (SQRT(DOT_PRODUCT(s, s)) < step_tol) THEN
                PRINT *, "  BFGS converged: step size too small"
                converged = .TRUE.
                DO i = 1, num_vars
                    x(i) = x_new(i)
                END DO
                f = f_new
                EXIT
            END IF
            
            s_dot_y = DOT_PRODUCT(s, y)
            
            IF (s_dot_y > 1e-10) THEN
                rho = 1.0 / s_dot_y
                
                ! FIX 4: Proper BFGS update using contiguous temporaries
                ! b_inv_new = (I - rho*s*y') * b_inv * (I - rho*y*s') + rho*s*s'
                
                ! First: temp_vec1 = b_inv * y
                CALL sgemv('N', num_vars, num_vars, 1.0, b_inv, num_vars, y, 1, 0.0, temp_vec1, 1)
                
                ! Compute y' * b_inv * y for later use
                y_dot_b_inv_y = DOT_PRODUCT(y, temp_vec1)
                
                ! Initialize b_inv_new = b_inv
                DO i = 1, num_vars
                    DO j = 1, num_vars
                        b_inv_new(i,j) = b_inv(i,j)
                    END DO
                END DO
                
                ! b_inv_new -= rho * s * (b_inv * y)'
                CALL sger(num_vars, num_vars, -rho, s, 1, temp_vec1, 1, b_inv_new, num_vars)
                
                ! temp_vec2 = b_inv' * y (for symmetric matrix, b_inv' * y = b_inv * y = temp_vec1)
                ! b_inv_new -= rho * (b_inv * y) * s'
                CALL sger(num_vars, num_vars, -rho, temp_vec1, 1, s, 1, b_inv_new, num_vars)
                
                ! b_inv_new += rho * rho * (y' * b_inv * y) * s * s'
                CALL sger(num_vars, num_vars, rho*rho*y_dot_b_inv_y, s, 1, s, 1, b_inv_new, num_vars)
                
                ! b_inv_new += rho * s * s'
                CALL sger(num_vars, num_vars, rho, s, 1, s, 1, b_inv_new, num_vars)
                
                ! Copy back
                DO i = 1, num_vars
                    DO j = 1, num_vars
                        b_inv(i,j) = b_inv_new(i,j)
                    END DO
                END DO
            ELSE
                PRINT *, "  Warning: Curvature condition violated, skipping BFGS update"
            END IF
            
            DO i = 1, num_vars
                x(i) = x_new(i)
                grad(i) = grad_new(i)
            END DO
            f = f_new
            
            IF (MOD(iter, 10) == 0) THEN
                PRINT *, "  BFGS iter", iter, ": fitness =", f, ", |grad| =", grad_norm, &
                        ", alpha =", alpha
            END IF
        END DO
        
        IF (.NOT. converged .AND. iter >= max_iter) THEN
            PRINT *, "  BFGS reached maximum iterations"
        END IF
        
        ! Explicit copy back
        DO i = 1, num_vars
            x_start(i) = x(i)
        END DO
        fitness_val = f
        
        PRINT *, "  BFGS finished: fitness =", fitness_val, " at x =", x_start
        PRINT *, "  Improvement:", f - f_origin
        PRINT *, "  Distance moved:", SQRT(SUM((x - x_origin)**2))
        PRINT *
        
    END SUBROUTINE refine_niche_bfgs
    
    SUBROUTINE compute_gradient_fd(x, grad, num_vars, fitness_func_id, lower, upper, num_eval)
        USE fitness_functions
        REAL, DIMENSION(num_vars), INTENT(IN) :: x
        REAL, DIMENSION(num_vars), INTENT(OUT) :: grad
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        REAL, INTENT(IN) :: lower, upper
        INTEGER, INTENT(INOUT) :: num_eval
        
        REAL, PARAMETER :: h = 1e-5
        REAL :: f_plus, f_minus
        REAL, DIMENSION(num_vars) :: x_temp
        INTEGER :: i, j
        
        DO i = 1, num_vars
            ! Explicit copy
            DO j = 1, num_vars
                x_temp(j) = x(j)
            END DO
            
            x_temp(i) = MIN(upper, x(i) + h)
            f_plus = evaluate_fitness_single(x_temp, num_vars, fitness_func_id)
            num_eval = num_eval + 1
            
            x_temp(i) = MAX(lower, x(i) - h)
            f_minus = evaluate_fitness_single(x_temp, num_vars, fitness_func_id)
            num_eval = num_eval + 1
            
            grad(i) = (f_plus - f_minus) / (2.0 * h)
        END DO
        
    END SUBROUTINE compute_gradient_fd
    
    SUBROUTINE line_search_backtracking(x, search_dir, f, grad, alpha, x_new, f_new, &
                                       num_vars, fitness_func_id, lower, upper, num_eval)
        USE fitness_functions
        REAL, DIMENSION(num_vars), INTENT(IN) :: x, search_dir, grad
        REAL, INTENT(IN) :: f
        REAL, INTENT(OUT) :: alpha, f_new
        REAL, DIMENSION(num_vars), INTENT(OUT) :: x_new
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        REAL, INTENT(IN) :: lower, upper
        INTEGER, INTENT(INOUT) :: num_eval
        
        REAL, PARAMETER :: c1 = 1e-4
        REAL, PARAMETER :: rho = 0.5
        REAL, PARAMETER :: alpha_init = 1.0
        INTEGER, PARAMETER :: max_iter = 20
        
        REAL :: directional_deriv
        INTEGER :: iter, i
        
        alpha = alpha_init
        directional_deriv = DOT_PRODUCT(grad, search_dir)
        
        DO iter = 1, max_iter
            DO i = 1, num_vars
                x_new(i) = x(i) + alpha * search_dir(i)
                x_new(i) = MAX(lower, MIN(upper, x_new(i)))
            END DO
            
            f_new = evaluate_fitness_single(x_new, num_vars, fitness_func_id)
            num_eval = num_eval + 1
            
            IF (f_new >= f + c1 * alpha * directional_deriv) THEN
                RETURN
            END IF
            
            alpha = rho * alpha
        END DO
        
        alpha = 1e-6
        DO i = 1, num_vars
            x_new(i) = x(i) + alpha * search_dir(i)
            x_new(i) = MAX(lower, MIN(upper, x_new(i)))
        END DO
        f_new = evaluate_fitness_single(x_new, num_vars, fitness_func_id)
        num_eval = num_eval + 1
        
    END SUBROUTINE line_search_backtracking
    
    FUNCTION evaluate_fitness_single(x, num_vars, fitness_func_id) RESULT(f)
        USE fitness_functions
        REAL, DIMENSION(num_vars), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        REAL :: f
        
        SELECT CASE(fitness_func_id)
            CASE(1)
                f = fitness_five_peaks(x, num_vars)
            CASE(2)
                f = fitness_rastrigin(x, num_vars)
            CASE(3)
                f = fitness_sphere(x, num_vars)
            CASE(4)
                f = fitness_schwefel(x, num_vars)
            CASE(5)
                f = fitness_equal_peaks(x, num_vars)
            CASE DEFAULT
                f = fitness_five_peaks(x, num_vars)
        END SELECT
        
    END FUNCTION evaluate_fitness_single
    
    SUBROUTINE consolidate_converged_niches(centers, center_fits, radius_adaptive, &
                                           num_nich, dims, merge_threshold_factor)
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        REAL, DIMENSION(num_nich), INTENT(INOUT) :: center_fits, radius_adaptive
        INTEGER, INTENT(IN) :: num_nich, dims
        REAL, INTENT(IN) :: merge_threshold_factor
        
        INTEGER :: i, j, k, merged_count, survivor_idx, merged_idx
        REAL :: dist, fitness_diff, merge_threshold
        LOGICAL :: should_merge
        REAL, PARAMETER :: fitness_tol = 0.01
        REAL, PARAMETER :: invalid_niche = -999.0
        
        merged_count = 0
        
        PRINT *, "Consolidating niches that converged to same optimum..."
        PRINT *
        
        DO i = 1, num_nich - 1
            IF (center_fits(i) < -900.0) CYCLE
            
            DO j = i + 1, num_nich
                IF (center_fits(j) < -900.0) CYCLE
                
                ! Explicit distance calculation
                dist = 0.0
                DO k = 1, dims
                    dist = dist + (centers(i,k) - centers(j,k))**2
                END DO
                dist = SQRT(dist)
                
                fitness_diff = ABS(center_fits(i) - center_fits(j))
                
                merge_threshold = merge_threshold_factor * MAX(radius_adaptive(i), radius_adaptive(j))
                
                should_merge = (dist < merge_threshold) .AND. &
                              (fitness_diff < fitness_tol * MAX(ABS(center_fits(i)), ABS(center_fits(j))))
                
                IF (should_merge) THEN
                    IF (center_fits(i) >= center_fits(j)) THEN
                        survivor_idx = i
                        merged_idx = j
                    ELSE
                        survivor_idx = j
                        merged_idx = i
                    END IF
                    
                    PRINT *, "  Merging duplicate niches", i, "and", j
                    PRINT *, "    Niche", i, ": pos =", (centers(i,k), k=1,dims), ", fitness =", center_fits(i)
                    PRINT *, "    Niche", j, ": pos =", (centers(j,k), k=1,dims), ", fitness =", center_fits(j)
                    PRINT *, "    Distance =", dist, ", threshold =", merge_threshold
                    PRINT *, "    Keeping niche", survivor_idx
                    PRINT *
                    
                    ! Invalidate merged niche
                    DO k = 1, dims
                        centers(merged_idx, k) = invalid_niche
                    END DO
                    center_fits(merged_idx) = invalid_niche
                    
                    merged_count = merged_count + 1
                    
                    IF (survivor_idx /= i) EXIT
                END IF
            END DO
        END DO
        
        IF (merged_count > 0) THEN
            PRINT *, "Consolidated", merged_count, " duplicate niche(s)"
        ELSE
            PRINT *, "No duplicate niches found - all converged to distinct optima"
        END IF
        PRINT *
        
    END SUBROUTINE consolidate_converged_niches

END MODULE bfgs_refinement

! ==============================================================================

PROGRAM vanesa
    USE blas_external
    USE fitness_functions
    USE bfgs_refinement
    IMPLICIT NONE
    
    REAL, PARAMETER :: invalid_niche = -999.0
    REAL, PARAMETER :: invalid_check = -998.0
    
    ! === CONFIGURABLE PARAMETERS (load from file in production) ===
    INTEGER, PARAMETER :: pop_size = 500
    INTEGER, PARAMETER :: max_gen = 30000  
    INTEGER, PARAMETER :: num_vars = 2 
    INTEGER, PARAMETER :: fitness_function = 1
    REAL, PARAMETER :: min_bound = -40.0
    REAL, PARAMETER :: max_bound = 40.0
    
    INTEGER, PARAMETER :: num_niches = 10
    INTEGER, PARAMETER :: min_niche_size = 4
    REAL, PARAMETER :: radius_factor = 0.10
    REAL, PARAMETER :: merge_radius_factor = 0.75
    INTEGER, PARAMETER :: merge_interval = 2
    REAL, PARAMETER :: consolidate_threshold = 0.2
    
    REAL, PARAMETER :: cmaes_c_mu = 0.5
    REAL, PARAMETER :: cmaes_damp = 2.0
    
    REAL, PARAMETER :: diversity_rate = 0.10
    REAL, PARAMETER :: reg_param = 1e-6
    INTEGER, PARAMETER :: print_interval = 5
    
    ! === STATE VARIABLES ===
    REAL, DIMENSION(num_niches) :: cmaes_sigma
    REAL, DIMENSION(num_niches, num_vars, num_vars) :: cmaes_cov
    REAL, DIMENSION(num_niches, num_vars) :: cmaes_pc, cmaes_ps
    REAL, DIMENSION(num_niches) :: cmaes_chin
    LOGICAL, DIMENSION(num_niches) :: cmaes_initialized
    
    REAL, DIMENSION(pop_size, num_vars) :: population, new_population
    REAL, DIMENSION(pop_size) :: fitness
    INTEGER, DIMENSION(pop_size) :: niche_membership
    
    REAL :: chin_calc
    REAL, DIMENSION(num_niches, num_vars) :: niche_centers
    REAL, DIMENSION(num_niches) :: niche_fitness
    REAL, DIMENSION(num_niches, num_vars, num_vars) :: niche_inv_cov
    LOGICAL, DIMENSION(num_niches) :: cov_valid
    REAL, DIMENSION(num_niches) :: niche_radius_adaptive
    REAL :: niche_radius_base
    
    INTEGER :: i, gen, num_eval, niche_id, member_count
    INTEGER, DIMENSION(pop_size) :: niche_members
    REAL :: time_start, time_end
    
    ! FIX 5: Proper sized work matrices
    REAL, DIMENSION(num_vars, num_vars) :: work_mat
    REAL, DIMENSION(num_vars) :: work_vec
    REAL :: trust_region
    
    ! === INITIALIZATION ===
    CALL CPU_TIME(time_start)
    CALL RANDOM_SEED()
    
    niche_radius_base = radius_factor * (max_bound - min_bound)
    niche_radius_adaptive = niche_radius_base
    
    CALL initialize_population_lhs(population, pop_size, num_vars, min_bound, max_bound)
    
    DO i = 1, num_niches
        CALL reset_cmaes_state(i, cmaes_sigma, cmaes_cov, cmaes_pc, cmaes_ps, &
                              cmaes_initialized, num_niches, num_vars)
    END DO
    
    cmaes_chin = SQRT(chi_square_critical(num_vars, 0.5))
    chin_calc = chi_square_critical(num_vars, 0.95)
    
    niche_centers = invalid_niche
    niche_fitness = invalid_niche
    cov_valid = .FALSE.
    num_eval = 0
    
    ! === MAIN EVOLUTIONARY LOOP ===
    DO gen = 1, max_gen
        print*, "1"
        CALL evaluate_fitness(population, fitness, pop_size, num_vars, num_eval)
        
        print*, "2"
        CALL update_niches_improved(population, fitness, niche_centers, niche_fitness, &
                                    num_niches, pop_size, num_vars, niche_radius_base)
        
        print*, "3"
        IF (gen == 1) THEN
            CALL assign_niches_euclidean(population, niche_centers, niche_membership, &
                                        pop_size, num_niches, num_vars, niche_radius_adaptive)
        END IF
        
        print*, "4"
        CALL update_niche_covariances(population, niche_membership, niche_centers, &
                                      niche_inv_cov, cov_valid, &
                                      pop_size, num_niches, num_vars)
        
        print*, "5"
        CALL assign_niches_mahalanobis(population, niche_centers, niche_inv_cov, &
                                      cov_valid, niche_membership, pop_size, &
                                      num_niches, num_vars, chin_calc)
        
        print*, "6"
        CALL update_adaptive_niche_radius(population, niche_membership, niche_centers, &
                                         niche_radius_adaptive, niche_radius_base, &
                                         pop_size, num_niches, num_vars)
        
        print*, "7"
        CALL create_new_generation_cmaes(new_population, niche_membership, &
                                        niche_centers, pop_size, num_vars, num_niches, &
                                        cmaes_sigma, cmaes_cov, min_bound, max_bound)
        
        print*, "8"
        CALL inject_diversity(new_population, niche_membership, pop_size, num_vars, &
                             min_bound, max_bound, diversity_rate)
        
        population = new_population
        
        print*, "9"
        DO niche_id = 1, num_niches
            IF (is_niche_valid(niche_id, niche_centers, num_niches, num_vars)) THEN
                CALL get_niche_members(niche_membership, niche_id, niche_members, &
                                      member_count, pop_size)
                
                IF (member_count >= min_niche_size) THEN
                    CALL update_cmaes_state_niche(population, fitness, &
                                                 niche_members, member_count, &
                                                 niche_centers, niche_id, &
                                                 cmaes_sigma(niche_id), &
                                                 cmaes_cov, niche_id, &
                                                 cmaes_pc, niche_id, &
                                                 cmaes_ps, niche_id, &
                                                 cmaes_chin(niche_id), &
                                                 num_vars, num_niches, pop_size, &
                                                 cmaes_initialized(niche_id))
                END IF
            END IF
        END DO
        
        print*, "10"
        IF (MOD(gen, merge_interval) == 0) THEN
            CALL merge_overlapping_niches(niche_centers, niche_fitness, &
                                         niche_radius_adaptive, cmaes_sigma, &
                                         cmaes_cov, cmaes_pc, cmaes_ps, cmaes_initialized, &
                                         num_niches, num_vars)
            CALL fix_ghost_niches(niche_centers, niche_membership, niche_radius_adaptive, &
                                 cmaes_sigma, pop_size, num_niches, num_vars)
        END IF
        
        print*, "11"
        IF (MOD(gen, print_interval) == 0) THEN
            CALL print_progress(gen, niche_centers, niche_fitness, niche_radius_adaptive, &
                               cmaes_sigma, num_eval, num_niches, num_vars)
        END IF
        print*, "12"
    END DO
    
    ! === BFGS REFINEMENT PHASE ===
    PRINT *, "==================================================="
    PRINT *, "Starting BFGS refinement phase"
    PRINT *, "==================================================="
    PRINT *
    
    DO niche_id = 1, num_niches
        IF (is_niche_valid(niche_id, niche_centers, num_niches, num_vars)) THEN
            PRINT *, "Refining niche", niche_id
            
            ! FIX 6: Explicit copy to avoid array slice issues
            DO i = 1, num_vars
                work_vec(i) = niche_centers(niche_id, i)
            END DO
            
            ! FIX 7: Explicit copy of covariance matrix slice
            DO j = 1, num_vars
                DO i = 1, num_vars
                    work_mat(i, j) = cmaes_sigma(niche_id)**2 * cmaes_cov(niche_id, i, j)
                END DO
            END DO
            
            trust_region = 1.5 * niche_radius_adaptive(niche_id)
            
            CALL refine_niche_bfgs(work_vec, work_mat, &
                                   niche_fitness(niche_id), num_vars, &
                                   fitness_function, min_bound, max_bound, &
                                   trust_region, num_eval)
            
            ! Copy back
            DO i = 1, num_vars
                niche_centers(niche_id, i) = work_vec(i)
            END DO
        END IF
    END DO
    
    PRINT *, "==================================================="
    PRINT *, "BFGS refinement completed"
    PRINT *, "==================================================="
    PRINT *
    
    CALL consolidate_converged_niches(niche_centers, niche_fitness, &
                                     niche_radius_adaptive, &
                                     num_niches, num_vars, consolidate_threshold)
    
    ! === FINALIZATION ===
    CALL print_final_results(niche_centers, niche_fitness, niche_radius_adaptive, &
                            cmaes_sigma, num_niches, num_vars)
    
    CALL CPU_TIME(time_end)
    PRINT *, "Time taken: ", time_end - time_start, " seconds"
    PRINT *, "Total evaluations: ", num_eval    
    
CONTAINS

    ! ==========================================================================
    ! UTILITY FUNCTIONS
    ! ==========================================================================
    
    FUNCTION is_niche_valid(niche_id, centers, num_nich, dims) RESULT(valid)
        INTEGER, INTENT(IN) :: niche_id, num_nich, dims
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        LOGICAL :: valid
        ! Better looking function :-)
        valid = .FALSE.
        IF (niche_id .gt. 0) valid = (centers(niche_id, 1) > invalid_check)

    END FUNCTION is_niche_valid
    
    FUNCTION euclidean_distance_explicit(point1, point2, dims) RESULT(dist)
        INTEGER, INTENT(IN) :: dims
        REAL, DIMENSION(dims), INTENT(IN) :: point1, point2
        REAL :: dist
        INTEGER :: i
        dist = 0.0
        DO i = 1, dims
            dist = dist + (point1(i) - point2(i))**2
        END DO
        dist = SQRT(dist)
    END FUNCTION euclidean_distance_explicit
    
    FUNCTION count_valid_niches(center_fits, num_nich) RESULT(count)
        REAL, DIMENSION(num_nich), INTENT(IN) :: center_fits
        INTEGER, INTENT(IN) :: num_nich
        INTEGER :: count, i
        count = 0
        DO i = 1, num_nich
            IF (center_fits(i) > invalid_niche) count = count + 1
        END DO
    END FUNCTION count_valid_niches
    
    SUBROUTINE swap_real(a, b)
        REAL, INTENT(INOUT) :: a, b
        REAL :: temp
        temp = a
        a = b
        b = temp
    END SUBROUTINE swap_real
    
    SUBROUTINE swap_int(a, b)
        INTEGER, INTENT(INOUT) :: a, b
        INTEGER :: temp
        temp = a
        a = b
        b = temp
    END SUBROUTINE swap_int
    
    FUNCTION inverse_normal_cdf(p) RESULT(z)
        REAL, INTENT(IN) :: p
        REAL :: z
        REAL :: q, r
        REAL, PARAMETER :: p_low  = 0.02425
        REAL, PARAMETER :: p_high = 0.97575
        
        REAL, PARAMETER :: a(6) = (/ -3.969683028665376e+01,  2.209460984245205e+02, &
                                     -2.759285104469687e+02,  1.383577518672690e+02, &
                                     -3.066479806614716e+01,  2.506628277459239e+00 /)
        REAL, PARAMETER :: b(5) = (/ -5.447609879822406e+01,  1.615858368580409e+02, &
                                     -1.556989798598866e+02,  6.680131188771972e+01, &
                                     -1.328068155288572e+01 /)
        
        REAL, PARAMETER :: c(6) = (/ -7.784894002430293e-03, -3.223964580411365e-01, &
                                     -2.400758277161838e+00, -2.549732539343734e+00, &
                                      4.374664141464968e+00,  2.938163982698783e+00 /)
        REAL, PARAMETER :: d(4) = (/  7.784695709041462e-03,  3.224671290700398e-01, &
                                      2.445134137142996e+00,  3.754408661907416e+00 /)
        
        IF (p < p_low) THEN
            q = SQRT(-2.0 * LOG(p))
            z = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
                ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0)
        ELSE IF (p <= p_high) THEN
            q = p - 0.5
            r = q * q
            z = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6)) * q / &
                (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0)
        ELSE
            q = SQRT(-2.0 * LOG(1.0 - p))
            z = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
                 ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0)
        END IF
        
    END FUNCTION inverse_normal_cdf
    
    FUNCTION chi_square_critical(df, percentile) RESULT(chi_crit)
        INTEGER, INTENT(IN) :: df
        REAL, INTENT(IN) :: percentile
        REAL :: chi_crit
        REAL :: z_alpha
        
        z_alpha = inverse_normal_cdf(percentile)
        chi_crit = df * (1.0 - 2.0/(9.0*df) + z_alpha*SQRT(2.0/(9.0*df)))**3
    END FUNCTION chi_square_critical

    ! ==========================================================================
    ! NICHE MANAGEMENT
    ! ==========================================================================

    SUBROUTINE get_niche_members(membership, niche_id, members, count, size)
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        INTEGER, INTENT(IN) :: niche_id, size
        INTEGER, DIMENSION(size), INTENT(OUT) :: members
        INTEGER, INTENT(OUT) :: count
        INTEGER :: i
        
        count = 0
        DO i = 1, size
            IF (membership(i) == niche_id) THEN
                count = count + 1
                members(count) = i
            END IF
        END DO
    END SUBROUTINE get_niche_members
    
    SUBROUTINE update_niches_improved(pop, fit, centers, center_fits, num_nich, size, dims, radius)
        REAL, DIMENSION(size, dims), INTENT(IN) :: pop
        REAL, DIMENSION(size), INTENT(IN) :: fit
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        REAL, DIMENSION(num_nich), INTENT(INOUT) :: center_fits
        INTEGER, INTENT(IN) :: num_nich, size, dims
        REAL, INTENT(IN) :: radius
        
        INTEGER :: i, j, k, idx, valid_niches
        INTEGER, DIMENSION(size) :: sorted_indices
        REAL :: dist, min_niche_separation
        LOGICAL :: is_new_niche
        REAL, DIMENSION(dims) :: temp_point1, temp_point2
        
        CALL sort_by_fitness(fit, sorted_indices, size)
        valid_niches = count_valid_niches(center_fits, num_nich)
        min_niche_separation = radius * 1.00
        
        DO i = 1, size
            idx = sorted_indices(i)
            is_new_niche = .TRUE.
            
            ! Explicit copy for distance calculation
            DO k = 1, dims
                temp_point1(k) = pop(idx, k)
            END DO
            
            DO j = 1, num_nich
                IF (.NOT. is_niche_valid(j, centers, num_nich, dims)) CYCLE
                
                DO k = 1, dims
                    temp_point2(k) = centers(j, k)
                END DO
                
                dist = euclidean_distance_explicit(temp_point1, temp_point2, dims)
                
                IF (dist < min_niche_separation) THEN
                    is_new_niche = .FALSE.
                    IF (fit(idx) > center_fits(j)) THEN
                        DO k = 1, dims
                            centers(j, k) = pop(idx, k)
                        END DO
                        center_fits(j) = fit(idx)
                    END IF
                    EXIT
                END IF
            END DO
            
            IF (is_new_niche .AND. valid_niches < num_nich) THEN
                DO j = 1, num_nich
                    IF (center_fits(j) < invalid_check) THEN
                        DO k = 1, dims
                            centers(j, k) = pop(idx, k)
                        END DO
                        center_fits(j) = fit(idx)
                        valid_niches = valid_niches + 1
                        EXIT
                    END IF
                END DO
            END IF
        END DO
    END SUBROUTINE update_niches_improved

    SUBROUTINE assign_niches_euclidean(pop, centers, membership, size, num_nich, dims, radius_adaptive)
        REAL, DIMENSION(size, dims), INTENT(IN) :: pop
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        INTEGER, DIMENSION(size), INTENT(OUT) :: membership
        INTEGER, INTENT(IN) :: size, num_nich, dims
        REAL, DIMENSION(num_nich), INTENT(IN) :: radius_adaptive
        INTEGER :: i, j, k, closest_niche
        REAL :: dist, min_dist
        REAL, DIMENSION(dims) :: temp_point1, temp_point2
        
        DO i = 1, size
            membership(i) = 0
            min_dist = HUGE(min_dist)
            closest_niche = 0
            
            DO k = 1, dims
                temp_point1(k) = pop(i, k)
            END DO
            
            DO j = 1, num_nich
                IF (.NOT. is_niche_valid(j, centers, num_nich, dims)) CYCLE
                
                DO k = 1, dims
                    temp_point2(k) = centers(j, k)
                END DO
                
                dist = euclidean_distance_explicit(temp_point1, temp_point2, dims) / radius_adaptive(j)
                
                IF (dist < min_dist) THEN
                    min_dist = dist
                    closest_niche = j
                END IF
            END DO
            
            IF (closest_niche > 0 .AND. min_dist <= 2.0) THEN
                membership(i) = closest_niche
            END IF
        END DO
    END SUBROUTINE assign_niches_euclidean

    SUBROUTINE assign_niches_mahalanobis(pop, centers, inv_cov, cov_valid, membership, &
                                        size, num_nich, dims, chi_threshold)
        REAL, DIMENSION(size, dims), INTENT(IN) :: pop
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        REAL, DIMENSION(num_nich, dims, dims), INTENT(IN) :: inv_cov
        LOGICAL, DIMENSION(num_nich), INTENT(IN) :: cov_valid
        INTEGER, DIMENSION(size), INTENT(OUT) :: membership
        INTEGER, INTENT(IN) :: size, num_nich, dims
        REAL, INTENT(IN) :: chi_threshold
        INTEGER :: i, j, k, closest_niche
        REAL :: dist, min_dist
        REAL, DIMENSION(dims) :: temp_point1, temp_point2
        REAL, DIMENSION(dims,dims) :: temp_inv_cov
        
        DO i = 1, size
            membership(i) = 0
            min_dist = HUGE(min_dist)
            closest_niche = 0
            
            DO k = 1, dims
                temp_point1(k) = pop(i, k)
            END DO
            
            DO j = 1, num_nich
                IF (.NOT. is_niche_valid(j, centers, num_nich, dims)) CYCLE
                
                DO k = 1, dims
                    temp_point2(k) = centers(j, k)
                END DO
                
                IF (cov_valid(j)) THEN
                    ! FIX 8: Explicit copy to avoid array slice issues
                    DO k = 1, dims
                        temp_inv_cov(k,1:dims) = inv_cov(j,k,1:dims)
                    END DO
                    dist = mahalanobis_distance_explicit(temp_point1, temp_point2, &
                                                        temp_inv_cov, dims)
                    dist = SQRT(2.0) * dist / SQRT(chi_threshold)
                ELSE
                    dist = euclidean_distance_explicit(temp_point1, temp_point2, dims)
                END IF
          
                IF (dist < min_dist) THEN
                    min_dist = dist
                    closest_niche = j
                END IF
            END DO
            
            IF (closest_niche > 0 .AND. min_dist <= 2.0) THEN
                membership(i) = closest_niche
            END IF
        END DO
    END SUBROUTINE assign_niches_mahalanobis
    
    SUBROUTINE update_adaptive_niche_radius(pop, membership, centers, radius_adaptive, radius_base, &
                                           size, num_nich, dims)
        REAL, DIMENSION(size, dims), INTENT(IN) :: pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        REAL, DIMENSION(num_nich), INTENT(OUT) :: radius_adaptive
        REAL, INTENT(IN) :: radius_base
        INTEGER, INTENT(IN) :: size, num_nich, dims
        
        INTEGER :: niche_id, i, j, count
        REAL :: mean_dist, dist
        REAL, DIMENSION(dims) :: temp_point1, temp_point2
        
        DO niche_id = 1, num_nich
            IF (.NOT. is_niche_valid(niche_id, centers, num_nich, dims)) THEN
                radius_adaptive(niche_id) = radius_base
                CYCLE
            END IF

            count = 0
            mean_dist = 0.0
            
            DO j = 1, dims
                temp_point2(j) = centers(niche_id, j)
            END DO
            
            DO i = 1, size
                IF (membership(i) == niche_id) THEN
                    DO j = 1, dims
                        temp_point1(j) = pop(i, j)
                    END DO
                    dist = euclidean_distance_explicit(temp_point1, temp_point2, dims)
                    mean_dist = mean_dist + dist
                    count = count + 1
                END IF
            END DO
            
            IF (count > 0) THEN
                mean_dist = mean_dist / REAL(count)
                radius_adaptive(niche_id) = MAX(0.05 * radius_base, 1.5 * mean_dist)
            END IF
        END DO
    END SUBROUTINE update_adaptive_niche_radius
    
    SUBROUTINE merge_overlapping_niches(centers, center_fits, radius_adaptive, &
                                       sigma, c, pc, ps, initialized, num_nich, dims)
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        REAL, DIMENSION(num_nich), INTENT(INOUT) :: center_fits, sigma, radius_adaptive
        REAL, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: c
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: pc, ps
        LOGICAL, DIMENSION(num_nich), INTENT(INOUT) :: initialized
        INTEGER, INTENT(IN) :: num_nich, dims
        
        INTEGER :: i, j, k, merged_count, survivor_idx, merged_idx
        REAL :: euclidean_dist, avg_radius, merge_threshold, new_radius
        REAL, DIMENSION(dims) :: temp_point1, temp_point2
        
        merged_count = 0
        
        DO i = 1, num_nich - 1
            IF (.NOT. is_niche_valid(i, centers, num_nich, dims)) CYCLE
            
            DO j = i + 1, num_nich
                IF (.NOT. is_niche_valid(j, centers, num_nich, dims)) CYCLE
                
                DO k = 1, dims
                    temp_point1(k) = centers(i, k)
                    temp_point2(k) = centers(j, k)
                END DO
                
                euclidean_dist = euclidean_distance_explicit(temp_point1, temp_point2, dims)
                avg_radius = (radius_adaptive(i) + radius_adaptive(j)) / 2.0
                merge_threshold = merge_radius_factor * avg_radius
                
                IF (euclidean_dist < merge_threshold) THEN
                    IF (center_fits(i) >= center_fits(j)) THEN
                        survivor_idx = i
                        merged_idx = j
                    ELSE
                        survivor_idx = j
                        merged_idx = i
                    END IF
                    
                    new_radius = MAX(radius_adaptive(survivor_idx), &
                                    euclidean_dist + radius_adaptive(merged_idx))
                    radius_adaptive(survivor_idx) = new_radius
                   !I will put this behind an option later – it creates an unbearable amount of print-spam 
                   !PRINT *, "Merged ", i, j, "with centers:", (centers(i,k), k=1,dims), &
                   !        " and", (centers(j,k), k=1,dims)
                    
                    CALL invalidate_niche(merged_idx, centers, center_fits, sigma, c, pc, ps, &
                                        initialized, num_nich, dims)
                    
                    merged_count = merged_count + 1
                    IF (survivor_idx /= i) EXIT 
                END IF
            END DO
        END DO
        
        IF (merged_count > 0) THEN
            PRINT *, "Merged ", merged_count, " overlapping niches"
        END IF
        
    END SUBROUTINE merge_overlapping_niches
    
    SUBROUTINE fix_ghost_niches(centers, membership, radius_adaptive, sigma, &
                                size, num_nich, dims)
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        REAL, DIMENSION(num_nich), INTENT(IN) :: radius_adaptive
        REAL, DIMENSION(num_nich), INTENT(INOUT) :: sigma
        INTEGER, INTENT(IN) :: size, num_nich, dims
        
        INTEGER :: niche_id, member_count, i, k
        REAL :: new_sigma
        
        DO niche_id = 1, num_nich
            IF (.NOT. is_niche_valid(niche_id, centers, num_nich, dims)) CYCLE
            
            member_count = 0
            DO i = 1, size
                IF (membership(i) == niche_id) member_count = member_count + 1
            END DO
            
            IF (member_count == 0) THEN
                new_sigma = MAX(0.5, radius_adaptive(niche_id) * 0.5)
                PRINT *, "WARNING: Ghost niche ", niche_id, &
                         " at ", (centers(niche_id,k), k=1,dims), &
                         " has 0 members. Resetting sigma from ", sigma(niche_id), &
                         " to ", new_sigma
                sigma(niche_id) = new_sigma
            END IF
        END DO
    END SUBROUTINE fix_ghost_niches

    SUBROUTINE invalidate_niche(niche_idx, centers, center_fits, sigma, c, pc, ps, &
                               initialized, num_nich, dims)
        INTEGER, INTENT(IN) :: niche_idx, num_nich, dims
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        REAL, DIMENSION(num_nich), INTENT(INOUT) :: center_fits, sigma
        REAL, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: c
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: pc, ps
        LOGICAL, DIMENSION(num_nich), INTENT(INOUT) :: initialized
        INTEGER :: i
        
        DO i = 1, dims
            centers(niche_idx, i) = invalid_niche
        END DO
        center_fits(niche_idx) = invalid_niche
        CALL reset_cmaes_state(niche_idx, sigma, c, pc, ps, initialized, num_nich, dims)
    END SUBROUTINE invalidate_niche

    SUBROUTINE reset_cmaes_state(niche_id, sigma, c, pc, ps, initialized, num_nich, dims)
        INTEGER, INTENT(IN) :: niche_id, num_nich, dims
        REAL, DIMENSION(num_nich), INTENT(INOUT) :: sigma
        REAL, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: c
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: pc, ps
        LOGICAL, DIMENSION(num_nich), INTENT(INOUT) :: initialized
        INTEGER :: j, k
        
        sigma(niche_id) = (max_bound - min_bound) / 4.0
        DO j = 1, dims
            DO k = 1, dims
                c(niche_id, j, k) = 0.0
            END DO
            c(niche_id, j, j) = 1.0
            pc(niche_id, j) = 0.0
            ps(niche_id, j) = 0.0
        END DO
        initialized(niche_id) = .FALSE.
    END SUBROUTINE reset_cmaes_state

    ! ==========================================================================
    ! COVARIANCE AND DISTANCE
    ! ==========================================================================

    FUNCTION mahalanobis_distance_explicit(point1, point2, inv_cov, dims) RESULT(dist)
        INTEGER, INTENT(IN) :: dims
        REAL, DIMENSION(dims), INTENT(IN) :: point1, point2
        REAL, DIMENSION(dims, dims), INTENT(IN) :: inv_cov
        REAL :: dist
        REAL, DIMENSION(dims) :: diff, temp
        INTEGER :: i
      
        ! Ugly debug code
        DO i = 1, dims
            diff(i) = point1(i) - point2(i)
        END DO
        CALL sgemv('N', dims, dims, 1.0, inv_cov, dims, diff, 1, 0.0, temp, 1)
        dist = SQRT(MAX(0.0, DOT_PRODUCT(diff, temp)))
    END FUNCTION mahalanobis_distance_explicit

    SUBROUTINE update_niche_covariances(pop, membership, centers, &
                                      inv_cov_matrices, cov_valid, size, num_nich, dims)
        REAL, DIMENSION(size, dims), INTENT(IN) :: pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        REAL, DIMENSION(num_nich, dims, dims), INTENT(OUT) :: inv_cov_matrices
        LOGICAL, DIMENSION(num_nich), INTENT(OUT) :: cov_valid
        INTEGER, INTENT(IN) :: size, num_nich, dims
        
        INTEGER :: niche_id, i, j, k, member_count, info, idx
        REAL, DIMENSION(dims) :: mean_vec, diff
        REAL, DIMENSION(dims, dims) :: cov_mat, work_mat
        INTEGER, DIMENSION(size) :: members
        INTEGER, DIMENSION(dims) :: ipiv
        REAL, DIMENSION(4*dims) :: work
        REAL :: rcond
        INTEGER, DIMENSION(dims) :: iwork
        
        cov_valid = .FALSE.
        
        DO niche_id = 1, num_nich
            IF (.NOT. is_niche_valid(niche_id, centers, num_nich, dims)) CYCLE
            
            CALL get_niche_members(membership, niche_id, members, member_count, size)
            
            IF (member_count >= 3) THEN
                mean_vec = 0.0
                DO i = 1, member_count
                    idx = members(i)
                    DO j = 1, dims
                        mean_vec(j) = mean_vec(j) + pop(idx, j)
                    END DO
                END DO
                DO j = 1, dims
                    mean_vec(j) = mean_vec(j) / REAL(member_count)
                END DO
                
                cov_mat = 0.0
                DO i = 1, member_count
                    idx = members(i)
                    DO j = 1, dims
                        diff(j) = pop(idx, j) - mean_vec(j)
                    END DO
                    
                    DO j = 1, dims
                        DO k = 1, dims
                            cov_mat(j, k) = cov_mat(j, k) + diff(j) * diff(k)
                        END DO
                    END DO
                END DO
                
                DO j = 1, dims
                    DO k = 1, dims
                        cov_mat(j, k) = cov_mat(j, k) / REAL(member_count - 1)
                    END DO
                    cov_mat(j, j) = cov_mat(j, j) + reg_param
                END DO
                
                work_mat = cov_mat
                CALL sgetrf(dims, dims, work_mat, dims, ipiv, info)
                IF (info == 0) THEN
                    CALL sgetri(dims, work_mat, dims, ipiv, work, dims, info)
                    IF (info == 0) THEN
                        CALL sgecon('1', dims, work_mat, dims, 1.0, rcond, work, iwork, info)
                        IF (info == 0) THEN
                            DO j = 1, dims
                                DO k = 1, dims
                                    inv_cov_matrices(niche_id, j, k) = work_mat(j, k)
                                END DO
                            END DO
                            cov_valid(niche_id) = .TRUE.
                        END IF
                    END IF
                END IF
                
                IF (.NOT. cov_valid(niche_id)) THEN
                    work_mat = cov_mat
                    CALL spotrf('U', dims, work_mat, dims, info)
                    IF (info == 0) THEN
                        CALL spotri('U', dims, work_mat, dims, info)
                        IF (info == 0) THEN
                            DO j = 1, dims
                                DO k = j+1, dims
                                    work_mat(k, j) = work_mat(j, k)
                                END DO
                            END DO
                            DO j = 1, dims
                                DO k = 1, dims
                                    inv_cov_matrices(niche_id, j, k) = work_mat(j, k)
                                END DO
                            END DO
                            cov_valid(niche_id) = .TRUE.
                        END IF
                    END IF
                END IF
            END IF
        END DO
    END SUBROUTINE update_niche_covariances
    
    ! ==========================================================================
    ! CMA-ES OPERATIONS
    ! ==========================================================================

    SUBROUTINE update_cmaes_state_niche(population, fitness, members, &
                                       member_count, niche_centers, niche_id, sigma, &
                                       cmaes_cov, cov_niche_id, &
                                       cmaes_pc, pc_niche_id, &
                                       cmaes_ps, ps_niche_id, &
                                       chin, dims, num_nich, size, initialized)
        REAL, DIMENSION(size,dims), INTENT(IN) :: population
        REAL, DIMENSION(size), INTENT(IN) :: fitness
        INTEGER, DIMENSION(size), INTENT(IN) :: members
        INTEGER, INTENT(IN) :: member_count, dims, niche_id, cov_niche_id
        INTEGER, INTENT(IN) :: pc_niche_id, ps_niche_id, num_nich, size
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: niche_centers
        REAL, INTENT(INOUT) :: sigma
        REAL, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: cmaes_cov
        REAL, DIMENSION(num_nich, dims), INTENT(INOUT) :: cmaes_pc, cmaes_ps
        REAL, INTENT(IN) :: chin
        LOGICAL, INTENT(INOUT) :: initialized
        
        INTEGER :: i, j, k, idx
        REAL, DIMENSION(member_count) :: member_fitness
        REAL :: mu_eff, cs, cc, c1, cmu
        REAL, DIMENSION(dims) :: mean_step, z_mean, niche_center_copy
        ! FIX 9: Use proper contiguous temp arrays
        REAL, DIMENSION(dims) :: pc_temp, ps_temp
        REAL, DIMENSION(dims, dims) :: c_new
        REAL :: ps_norm, hsig
    
        
        DO i = 1, member_count
            member_fitness(i) = fitness(members(i))
        END DO
        
        IF (.NOT. initialized) THEN
            initialized = .TRUE.
            RETURN
        END IF
        
        mu_eff = REAL(MIN(member_count, dims)) / 2.0
        cs = (mu_eff + 2.0) / (dims + mu_eff + 5.0)
        cc = (4.0 + mu_eff/dims) / (dims + 4.0 + 2.0*mu_eff/dims)
        c1 = 2.0 / ((dims + 1.3)**2 + mu_eff)
        cmu = MIN(1.0 - c1, 2.0*(mu_eff - 2.0 + 1.0/mu_eff) / ((dims+2.0)**2 + mu_eff))
        
        ! Copy niche center
        DO i = 1, dims
            niche_center_copy(i) = niche_centers(niche_id, i)
        END DO
        
        mean_step = 0.0
        DO i = 1, MIN(INT(mu_eff), member_count)
            idx = members(i)
            DO j = 1, dims
                mean_step(j) = mean_step(j) + (population(idx, j) - niche_center_copy(j)) / sigma
            END DO
        END DO
        DO j = 1, dims
            mean_step(j) = mean_step(j) / REAL(MIN(INT(mu_eff), member_count))
        END DO
        
        ! Update ps - FIX 10: Copy to local contiguous array first
        DO i = 1, dims
            ps_temp(i) = cmaes_ps(ps_niche_id, i)
        END DO
        
        DO i = 1, dims
            ps_temp(i) = (1.0 - cs) * ps_temp(i) + SQRT(cs * (2.0 - cs)) * mean_step(i)
        END DO
        
        ps_norm = 0.0
        DO i = 1, dims
            ps_norm = ps_norm + ps_temp(i)**2
        END DO
        ps_norm = SQRT(ps_norm)
        
        hsig = 0.0
        IF (ps_norm / SQRT(1.0 - (1.0-cs)**2) < 1.4 * chin) THEN
            hsig = 1.0
        END IF
        
        ! Update pc - FIX 11: Copy to local contiguous array first
        DO i = 1, dims
            pc_temp(i) = cmaes_pc(pc_niche_id, i)
        END DO
        
        DO i = 1, dims
            pc_temp(i) = (1.0 - cc) * pc_temp(i) + hsig * SQRT(cc * (2.0 - cc)) * mean_step(i)
        END DO
        
        ! Update covariance
        DO k = 1, dims
            DO j = 1, dims
                c_new(j, k) = (1.0 - c1 - cmu) * cmaes_cov(cov_niche_id, j, k)
            END DO
        END DO
        
        ! FIX 12: Use local contiguous pc_temp for BLAS call
        CALL sger(dims, dims, c1, pc_temp, 1, pc_temp, 1, c_new, dims)
        
        DO i = 1, MIN(INT(mu_eff), member_count)
            idx = members(i)
            DO j = 1, dims
                z_mean(j) = (population(idx, j) - niche_center_copy(j)) / sigma
            END DO
            CALL sger(dims, dims, cmu/REAL(mu_eff), z_mean, 1, z_mean, 1, c_new, dims)
        END DO
        
        DO j = 1, dims
            DO k = 1, dims
                cmaes_cov(cov_niche_id, j, k) = c_new(j, k)
            END DO
        END DO
        
        DO i = 1, dims
            cmaes_cov(cov_niche_id, i, i) = cmaes_cov(cov_niche_id, i, i) + reg_param
        END DO
        
        ! Copy back pc and ps
        DO i = 1, dims
            cmaes_pc(pc_niche_id, i) = pc_temp(i)
            cmaes_ps(ps_niche_id, i) = ps_temp(i)
        END DO
        
        sigma = sigma * EXP((cs / cmaes_damp) * (ps_norm / chin - 1.0))
        sigma = MAX(1e-6 * (max_bound - min_bound), &
                   MIN(1.0 * (max_bound - min_bound), sigma))
        
    END SUBROUTINE update_cmaes_state_niche
    
    SUBROUTINE create_new_generation_cmaes(new_pop, membership, centers, &
                                           size, dims, num_nich, sigma_cmaes, c_cmaes, &
                                           lower, upper)
        REAL, DIMENSION(size, dims), INTENT(OUT) :: new_pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        INTEGER, INTENT(IN) :: size, dims, num_nich
        REAL, DIMENSION(num_nich), INTENT(IN) :: sigma_cmaes
        REAL, DIMENSION(num_nich, dims, dims), INTENT(IN) :: c_cmaes
        REAL, INTENT(IN) :: lower, upper
        
        INTEGER :: i, j, k, niche_id
        REAL, DIMENSION(dims) :: sample, center_copy
        REAL, DIMENSION(dims, dims) :: cov_copy
        
        DO i = 1, size
            niche_id = membership(i)
            
            IF (niche_id > 0 .AND. is_niche_valid(niche_id, centers, num_nich, dims)) THEN
                ! FIX 13: Explicit copy to ensure contiguous memory
                DO j = 1, dims
                    center_copy(j) = centers(niche_id, j)
                END DO
                
                DO j = 1, dims
                    DO k = 1, dims
                        cov_copy(j, k) = c_cmaes(niche_id, j, k)
                    END DO
                END DO
                
                CALL sample_from_cmaes_distribution(center_copy, &
                                                   sigma_cmaes(niche_id), &
                                                   cov_copy, sample, dims)
                DO j = 1, dims
                    new_pop(i, j) = MAX(lower, MIN(upper, sample(j)))
                END DO
            ELSE
                CALL random_uniform_vector(sample, dims, lower, upper)
                DO j = 1, dims
                    new_pop(i, j) = sample(j)
                END DO
            END IF
        END DO
        
    END SUBROUTINE create_new_generation_cmaes
    
    SUBROUTINE sample_from_cmaes_distribution(mean, sigma, c, sample, dims)
        REAL, DIMENSION(dims), INTENT(IN) :: mean
        REAL, INTENT(IN) :: sigma
        REAL, DIMENSION(dims, dims), INTENT(IN) :: c
        REAL, DIMENSION(dims), INTENT(OUT) :: sample
        INTEGER, INTENT(IN) :: dims
        
        REAL, DIMENSION(dims) :: z, y
        INTEGER :: i
        REAL :: u1, u2, radius, z1, z2
        REAL, PARAMETER :: two_pi = 6.283185307179586
        
        DO i = 1, dims, 2
            CALL RANDOM_NUMBER(u1)
            CALL RANDOM_NUMBER(u2)
            
            radius = SQRT(-2.0 * LOG(u1 + 1e-10))
            z1 = radius * COS(two_pi * u2)
            z2 = radius * SIN(two_pi * u2)
            
            z(i) = z1
            IF (i + 1 <= dims) THEN
                z(i + 1) = z2
            END IF
        END DO
        
        ! FIX 14: c is already a contiguous copy from caller
        CALL sgemv('N', dims, dims, 1.0, c, dims, z, 1, 0.0, y, 1)
        DO i = 1, dims
            sample(i) = mean(i) + sigma * y(i)
        END DO
    
    END SUBROUTINE sample_from_cmaes_distribution  

    ! ==========================================================================
    ! POPULATION MANAGEMENT
    ! ==========================================================================
    
    SUBROUTINE initialize_population_lhs(pop, size, dims, lower, upper)
        REAL, DIMENSION(size, dims), INTENT(OUT) :: pop
        INTEGER, INTENT(IN) :: size, dims
        REAL, INTENT(IN) :: lower, upper
        INTEGER :: i, j, k
        REAL :: r
        INTEGER, DIMENSION(size) :: perm
        
        DO j = 1, dims
            DO i = 1, size
                perm(i) = i
            END DO
            
            DO i = size, 2, -1
                CALL RANDOM_NUMBER(r)
                k = 1 + INT(r * i)
                CALL swap_int(perm(i), perm(k))
            END DO
            
            DO i = 1, size
                CALL RANDOM_NUMBER(r)
                pop(i, j) = lower + (upper - lower) * (perm(i) - 1 + r) / REAL(size)
            END DO
        END DO
    END SUBROUTINE initialize_population_lhs
    
    SUBROUTINE inject_diversity(pop, membership, size, dims, lower, upper, rate)
        REAL, DIMENSION(size, dims), INTENT(INOUT) :: pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        INTEGER, INTENT(IN) :: size, dims
        REAL, INTENT(IN) :: lower, upper, rate
        
        INTEGER :: i, j, num_inject
        REAL :: r
        REAL, DIMENSION(dims) :: temp_vec
        
        num_inject = INT(rate * REAL(size))
        
        DO i = 1, size
            IF (i <= num_inject) THEN
                CALL random_uniform_vector(temp_vec, dims, lower, upper)
                DO j = 1, dims
                    pop(i, j) = temp_vec(j)
                END DO
            ELSE IF (membership(i) == 0) THEN
                CALL RANDOM_NUMBER(r)
                IF (r < 0.3) THEN
                    CALL random_uniform_vector(temp_vec, dims, lower, upper)
                    DO j = 1, dims
                        pop(i, j) = temp_vec(j)
                    END DO
                END IF
            END IF
        END DO
    END SUBROUTINE inject_diversity
    
    SUBROUTINE random_uniform_vector(vec, dims, lower, upper)
        REAL, DIMENSION(dims), INTENT(OUT) :: vec
        INTEGER, INTENT(IN) :: dims
        REAL, INTENT(IN) :: lower, upper
        INTEGER :: i
        REAL :: r
        
        DO i = 1, dims
            CALL RANDOM_NUMBER(r)
            vec(i) = lower + r * (upper - lower)
        END DO
    END SUBROUTINE random_uniform_vector
    
    ! ==========================================================================
    ! FITNESS EVALUATION AND SORTING
    ! ==========================================================================
    
    SUBROUTINE evaluate_fitness(pop, fit, size, dims, evalfs)
        REAL, DIMENSION(size, dims), INTENT(IN) :: pop
        REAL, DIMENSION(size), INTENT(OUT) :: fit
        INTEGER, INTENT(IN) :: size, dims
        INTEGER, INTENT(INOUT) :: evalfs
        INTEGER :: i, j
        REAL, DIMENSION(dims) :: temp_point
        
        DO i = 1, size
            evalfs = evalfs + 1
            
            DO j = 1, dims
                temp_point(j) = pop(i, j)
            END DO
            
            SELECT CASE(fitness_function)
                CASE(1)
                    fit(i) = fitness_five_peaks(temp_point, dims)
                CASE(2)
                    fit(i) = fitness_rastrigin(temp_point, dims)
                CASE(3)
                    fit(i) = fitness_sphere(temp_point, dims)
                CASE(4)
                    fit(i) = fitness_schwefel(temp_point, dims)
                CASE(5)
                    fit(i) = fitness_equal_peaks(temp_point, dims)
                CASE DEFAULT
                    fit(i) = fitness_five_peaks(temp_point, dims)
            END SELECT
        END DO
    END SUBROUTINE evaluate_fitness
    
    SUBROUTINE sort_by_fitness(fitness, indices, size)
        INTEGER, INTENT(IN) :: size
        REAL, DIMENSION(size), INTENT(IN) :: fitness
        INTEGER, DIMENSION(size), INTENT(OUT) :: indices
        
        REAL, DIMENSION(size) :: fitness_copy
        INTEGER :: i
        
        DO i = 1, size
            indices(i) = i
            fitness_copy(i) = fitness(i)
        END DO
        
        CALL quicksort_descending(fitness_copy, indices, 1, size)
        
    END SUBROUTINE sort_by_fitness

    RECURSIVE SUBROUTINE quicksort_descending(fitness, indices, left, right)
        REAL, DIMENSION(:), INTENT(INOUT) :: fitness
        INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
        INTEGER, INTENT(IN) :: left, right
        
        INTEGER :: pivot_idx
        
        IF (left < right) THEN
            CALL partition_descending(fitness, indices, left, right, pivot_idx)
            CALL quicksort_descending(fitness, indices, left, pivot_idx - 1)
            CALL quicksort_descending(fitness, indices, pivot_idx + 1, right)
        END IF
        
    END SUBROUTINE quicksort_descending

    SUBROUTINE partition_descending(fitness, indices, left, right, pivot_idx)
        REAL, DIMENSION(:), INTENT(INOUT) :: fitness
        INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
        INTEGER, INTENT(IN) :: left, right
        INTEGER, INTENT(OUT) :: pivot_idx
        
        REAL :: pivot_value
        INTEGER :: i
        
        pivot_idx = (left + right) / 2
        pivot_value = fitness(pivot_idx)
        
        CALL swap_real(fitness(pivot_idx), fitness(right))
        CALL swap_int(indices(pivot_idx), indices(right))
        
        pivot_idx = left
        DO i = left, right - 1
            IF (fitness(i) > pivot_value .OR. &
                (ABS(fitness(i) - pivot_value) < 1e-10 .AND. indices(i) < indices(right))) THEN
                CALL swap_real(fitness(i), fitness(pivot_idx))
                CALL swap_int(indices(i), indices(pivot_idx))
                pivot_idx = pivot_idx + 1
            END IF
        END DO
        
        CALL swap_real(fitness(pivot_idx), fitness(right))
        CALL swap_int(indices(pivot_idx), indices(right))
        
    END SUBROUTINE partition_descending

    ! ==========================================================================
    ! OUTPUT
    ! ==========================================================================
    
    SUBROUTINE print_progress(gen, centers, center_fits, radius_adaptive, sigma, &
                             num_eval, num_nich, dims)
        INTEGER, INTENT(IN) :: gen, num_eval, num_nich, dims
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        REAL, DIMENSION(num_nich), INTENT(IN) :: center_fits, radius_adaptive, sigma
        INTEGER :: i, j
        
        PRINT *, "Generation ", gen
        PRINT *, "Best niches found:"
        DO i = 1, num_nich
            IF (center_fits(i) > invalid_niche) THEN
                PRINT *, "  Niche ", i, ": Position = ", (centers(i,j), j=1,dims), &
                     " Fitness = ", center_fits(i), &
                     " Peak width = ", radius_adaptive(i), &
                     " Sigma = ", sigma(i)
            END IF
        END DO
        PRINT *, "Total evaluations: ", num_eval
        PRINT *
    END SUBROUTINE print_progress
    
    SUBROUTINE print_final_results(centers, center_fits, radius_adaptive, sigma, num_nich, dims)
        REAL, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        REAL, DIMENSION(num_nich), INTENT(IN) :: center_fits, radius_adaptive, sigma
        INTEGER, INTENT(IN) :: num_nich, dims
        
        INTEGER :: i, j, k, valid_count
        INTEGER, DIMENSION(num_nich) :: sorted_indices
        REAL, DIMENSION(num_nich) :: temp_fitness
        
        PRINT *, "==================================================="
        PRINT *, "Optimization completed"
        PRINT *, "Final niches found (sorted by fitness):"
        PRINT *, "==================================================="
        
        valid_count = 0
        DO i = 1, num_nich
            IF (center_fits(i) > invalid_niche) THEN
                valid_count = valid_count + 1
                sorted_indices(valid_count) = i
                temp_fitness(valid_count) = center_fits(i)
            END IF
        END DO
        
        DO j = 1, valid_count - 1
            DO k = j + 1, valid_count
                IF (temp_fitness(j) < temp_fitness(k)) THEN
                    CALL swap_real(temp_fitness(j), temp_fitness(k))
                    CALL swap_int(sorted_indices(j), sorted_indices(k))
                END IF
            END DO
        END DO
        
        DO j = 1, valid_count
            i = sorted_indices(j)
            PRINT *, "Rank ", j, " - Niche ", i, ":"
            PRINT *, "  Position    = ", (centers(i,k), k=1,dims)
            PRINT *, "  Fitness     = ", center_fits(i)
            PRINT *, "  Peak width  = ", radius_adaptive(i)
            PRINT *, "  Final Sigma = ", sigma(i)
            PRINT *
        END DO
        
        PRINT *, "Total valid niches found: ", valid_count
        PRINT *, "==================================================="
    END SUBROUTINE print_final_results
    
END PROGRAM vanesa
