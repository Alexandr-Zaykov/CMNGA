! ============================================================== !
! VANESA: Variable Adaptive Niching Evolution Strategy Algorithm !
! ============================================================== !

MODULE blas_external
    IMPLICIT NONE
    EXTERNAL :: dgemv, dger, dgetrf, dgetri, dgecon, dpotrf, dpotri
END MODULE blas_external

MODULE fitness_functions
    IMPLICIT NONE
    
CONTAINS
    
    FUNCTION fitness_five_peaks(x, dims) RESULT(f)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION :: f
        DOUBLE PRECISION :: a(5), b(5), h(5), w(5)
        INTEGER :: j
        
        DATA a /-20.0D0, -5.0D0, 0.0D0, 30.0D0, 30.0D0/
        DATA b /-20.0D0, -25.0D0, 30.0D0, 0.0D0, -30.0D0/
        DATA h /0.4D0, 0.2D0, 0.7D0, 1.0D0, 0.05D0/
        DATA w /0.02D0, 0.5D0, 0.01D0, 2.0D0, 0.1D0/
        
        IF (dims == 1) THEN
            f = DSIN(5.0D0 * x(1))**6 * DEXP(-x(1)**2) + 1.0D0
        ELSE
            f = 0.0D0
            DO j = 1, 5
                f = f + h(j) / (1.0D0 + w(j) * ((x(1)-a(j))**2 + (x(2)-b(j))**2))
            END DO
            f = f + 1.0D0
            
            IF (dims > 2) THEN
                DO j = 3, dims
                    f = f * DEXP(-(x(j)**2)/4.0D0)
                END DO
            END IF
        END IF
    END FUNCTION fitness_five_peaks
    
    FUNCTION fitness_rastrigin(x, dims) RESULT(f)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION :: f
        INTEGER :: i
        DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846D0

        
        f = 10.0D0 * dims
        DO i = 1, dims
            f = f + x(i)**2 - 10.0D0 * DCOS(2.0D0 * pi * x(i))
        END DO
        ! Final f =~ 40.35*num_vars
        ! Maxima are at +-4.52 when using x e {–5..5}
       ! If I want to switch it to minimization.
       ! f = 400 - f
    END FUNCTION fitness_rastrigin

    FUNCTION fitness_slow_rastrigin(x, dims) RESULT(f)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION :: f, dummy_work
        INTEGER :: i,j
        DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846D0

        dummy_work = 0.0D0 
        f = 10.0D0 * dims
        DO i = 1, dims
            DO j = 1, 100  ! ← Control parameter (adjust from 10 to 1000)
                dummy_work = dummy_work + DSIN(x(i) * j * 0.01D0) * DCOS(x(i) * j * 0.01D0)
                dummy_work = dummy_work + DEXP(-DABS(x(i)) * 0.001D0 * j)
                dummy_work = dummy_work + DSQRT(DABS(x(i) * j * 0.01D0) + 1.0D0)
            END DO
        END DO
        DO i = 1, dims
            f = f + x(i)**2 - 10.0D0 * DCOS(2.0D0 * pi * x(i))
        END DO
        f = f + dummy_work*1.0D-12
        ! Final f =~ 40.35*num_vars
        ! Maxima are at +-4.52 when using x e {–5..5}
       ! If I want to switch it to minimization.
       ! f = 400 - f
    END FUNCTION fitness_slow_rastrigin
    
    FUNCTION fitness_sphere(x, dims) RESULT(f)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION :: f
        
        f = -SUM(x**2) + 100.0D0
    END FUNCTION fitness_sphere
    
    FUNCTION fitness_schwefel(x, dims) RESULT(f)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION :: f
        INTEGER :: i
        
        f = 0.0D0
        DO i = 1, dims
            f = f + x(i) * DSIN(DSQRT(DABS(x(i))))
        END DO
        f = 418.9829D0 * dims + f
    END FUNCTION fitness_schwefel
    
    FUNCTION fitness_equal_peaks(x, dims) RESULT(f)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION :: f
        DOUBLE PRECISION :: peaks(4, 2)
        INTEGER :: i
        DOUBLE PRECISION :: dist, peak_val
        
        DATA peaks /-20.0D0, -20.0D0, -20.0D0, 20.0D0, 20.0D0, -20.0D0, 20.0D0, 20.0D0/
        
        f = 0.0D0
        DO i = 1, 4
            dist = DSQRT((x(1) - peaks(i,1))**2 + (x(2) - peaks(i,2))**2)
            peak_val = DEXP(-dist**2 / 50.0D0)
            f = DMAX1(f, peak_val)
        END DO
        f = f + 1.0D0
    END FUNCTION fitness_equal_peaks
    
END MODULE fitness_functions

! ==============================================================================

MODULE bfgs_refinement
    USE blas_external
    IMPLICIT NONE
    
CONTAINS

    SUBROUTINE refine_niche_bfgs(x_start, b_inv_init, fitness_val, num_vars, &
                             fitness_func_id, lower, upper, num_eval, did_converge)
        DOUBLE PRECISION, DIMENSION(num_vars), INTENT(INOUT) :: x_start
        DOUBLE PRECISION, DIMENSION(num_vars, num_vars), INTENT(IN) :: b_inv_init
        DOUBLE PRECISION, INTENT(INOUT) :: fitness_val
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        DOUBLE PRECISION, INTENT(IN) :: lower, upper
        INTEGER, INTENT(INOUT) :: num_eval
        LOGICAL, INTENT(OUT) :: did_converge
        
        INTEGER, PARAMETER :: max_iter = 50
        DOUBLE PRECISION, PARAMETER :: grad_tol = 1.0D-5
        DOUBLE PRECISION, PARAMETER :: step_tol = 1.0D-8
        DOUBLE PRECISION, PARAMETER :: f_tol = 1.0D-8
        DOUBLE PRECISION, PARAMETER :: min_grad_for_refinement = 1.0D-4
        
        DOUBLE PRECISION, DIMENSION(num_vars) :: x, x_new, grad, grad_new, search_dir, s, y
        DOUBLE PRECISION, DIMENSION(num_vars) :: temp_vec1
        DOUBLE PRECISION, DIMENSION(num_vars, num_vars) :: b_inv, b_inv_new
        DOUBLE PRECISION :: f, f_new, alpha, grad_norm, s_dot_y, rho, y_dot_b_inv_y
        DOUBLE PRECISION :: improvement, f_origin
        INTEGER :: iter, i, j
        LOGICAL :: converged
        
        ! Explicit copy to avoid array slice issues
        DO i = 1, num_vars
            x(i) = x_start(i)
        END DO
        f = fitness_val
        f_origin = fitness_val
        
        ! Explicit copy of matrix to ensure contiguous memory
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
        
        grad_norm = DSQRT(DOT_PRODUCT(grad, grad))
        
        PRINT *, "  Starting BFGS refinement at x =", x
        PRINT *, "  Initial fitness =", f
        PRINT *, "  Initial gradient norm =", grad_norm
        
        IF (grad_norm < min_grad_for_refinement) THEN
            PRINT *, "  Skipping: Already converged (gradient too small)"
            PRINT *
            did_converge = .TRUE.
            RETURN
        END IF
        
        DO iter = 1, max_iter
            grad_norm = DSQRT(DOT_PRODUCT(grad, grad))
            
            IF (grad_norm < grad_tol) THEN
                PRINT *, "  BFGS converged: gradient tolerance reached"
                converged = .TRUE.
                EXIT
            END IF
            
            ! Use contiguous temp_vec1 for search direction
            CALL dgemv('N', num_vars, num_vars, -1.0D0, b_inv, num_vars, grad, 1, 0.0D0, search_dir, 1)
            
            CALL line_search_backtracking(x, search_dir, f, grad, alpha, x_new, f_new, &
                                         num_vars, fitness_func_id, lower, upper, num_eval)
            
                      
            improvement = f_new - f
            
            IF (improvement < -1.0D-6 .AND. iter > 3) THEN
                PRINT *, "  BFGS stopped: Fitness degrading (change =", improvement, ")"
                converged = .TRUE.
                EXIT
            END IF
            
            IF (DABS(improvement) < f_tol .AND. iter > 5) THEN
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
            
            IF (DSQRT(DOT_PRODUCT(s, s)) < step_tol) THEN
                PRINT *, "  BFGS converged: step size too small"
                converged = .TRUE.
                DO i = 1, num_vars
                    x(i) = x_new(i)
                END DO
                f = f_new
                EXIT
            END IF
            
            s_dot_y = DOT_PRODUCT(s, y)
            
            IF (s_dot_y > 1.0D-10) THEN
                rho = 1.0D0 / s_dot_y
                
                ! Proper BFGS update using contiguous temporaries
                ! b_inv_new = (I - rho*s*y') * b_inv * (I - rho*y*s') + rho*s*s'
                
                ! First: temp_vec1 = b_inv * y
                CALL dgemv('N', num_vars, num_vars, 1.0D0, b_inv, num_vars, y, 1, 0.0D0, temp_vec1, 1)
                
                ! Compute y' * b_inv * y for later use
                y_dot_b_inv_y = DOT_PRODUCT(y, temp_vec1)
                
                ! Initialize b_inv_new = b_inv
                DO i = 1, num_vars
                    DO j = 1, num_vars
                        b_inv_new(i,j) = b_inv(i,j)
                    END DO
                END DO
                
                ! b_inv_new -= rho * s * (b_inv * y)'
                CALL dger(num_vars, num_vars, -rho, s, 1, temp_vec1, 1, b_inv_new, num_vars)
                
                ! b_inv_new -= rho * (b_inv * y) * s'
                CALL dger(num_vars, num_vars, -rho, temp_vec1, 1, s, 1, b_inv_new, num_vars)
                
                ! b_inv_new += rho * rho * (y' * b_inv * y) * s * s'
                CALL dger(num_vars, num_vars, rho*rho*y_dot_b_inv_y, s, 1, s, 1, b_inv_new, num_vars)
                
                ! b_inv_new += rho * s * s'
                CALL dger(num_vars, num_vars, rho, s, 1, s, 1, b_inv_new, num_vars)
                
                ! Copy back
                DO j = 1, num_vars
                    DO i = 1, num_vars
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
            PRINT *, "  BFGS reached maximum iterations!"
        END IF
        
        ! Explicit copy back
        DO i = 1, num_vars
            x_start(i) = x(i)
        END DO
        fitness_val = f
        ! convergence flag
        did_converge = converged
        
        PRINT *, "  BFGS finished: fitness =", fitness_val, " at x =", x_start
        PRINT *, "  Improvement:", f - f_origin
        PRINT *
        
    END SUBROUTINE refine_niche_bfgs
    
    SUBROUTINE compute_gradient_fd(x, grad, num_vars, fitness_func_id, lower, upper, num_eval)
        USE fitness_functions
        DOUBLE PRECISION, DIMENSION(num_vars), INTENT(IN) :: x
        DOUBLE PRECISION, DIMENSION(num_vars), INTENT(OUT) :: grad
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        DOUBLE PRECISION, INTENT(IN) :: lower, upper
        INTEGER, INTENT(INOUT) :: num_eval
        
        DOUBLE PRECISION, PARAMETER :: h = 1.0D-5
        DOUBLE PRECISION :: f_plus, f_minus
        DOUBLE PRECISION, DIMENSION(num_vars) :: x_temp
        INTEGER :: i, j
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, x_temp, f_plus, f_minus) SCHEDULE(STATIC)
        DO i = 1, num_vars
            ! Explicit copy
            DO j = 1, num_vars
                x_temp(j) = x(j)
            END DO
            
            x_temp(i) = DMIN1(upper, x(i) + h)
            f_plus = evaluate_fitness_single(x_temp, num_vars, fitness_func_id)
            
            x_temp(i) = DMAX1(lower, x(i) - h)
            f_minus = evaluate_fitness_single(x_temp, num_vars, fitness_func_id)
            
            grad(i) = (f_plus - f_minus) / (2.0D0 * h)
        END DO
        !$OMP END PARALLEL DO
        
        ! Update evaluation counter after parallel region
        num_eval = num_eval + 2 * num_vars
        
    END SUBROUTINE compute_gradient_fd
    
    SUBROUTINE line_search_backtracking(x, search_dir, f, grad, alpha, x_new, f_new, &
                                       num_vars, fitness_func_id, lower, upper, num_eval)
        USE fitness_functions
        DOUBLE PRECISION, DIMENSION(num_vars), INTENT(IN) :: x, search_dir, grad
        DOUBLE PRECISION, INTENT(IN) :: f
        DOUBLE PRECISION, INTENT(OUT) :: alpha, f_new
        DOUBLE PRECISION, DIMENSION(num_vars), INTENT(OUT) :: x_new
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        DOUBLE PRECISION, INTENT(IN) :: lower, upper
        INTEGER, INTENT(INOUT) :: num_eval
        
        DOUBLE PRECISION, PARAMETER :: c1 = 1.0D-4
        DOUBLE PRECISION, PARAMETER :: rho = 0.5D0
        DOUBLE PRECISION, PARAMETER :: alpha_init = 1.0D0
        INTEGER, PARAMETER :: max_iter = 20
        
        DOUBLE PRECISION :: directional_deriv
        INTEGER :: iter, i
        
        alpha = alpha_init
        directional_deriv = DOT_PRODUCT(grad, search_dir)
        
        DO iter = 1, max_iter
            DO i = 1, num_vars
                x_new(i) = x(i) + alpha * search_dir(i)
                x_new(i) = DMAX1(lower, DMIN1(upper, x_new(i)))
            END DO
            
            f_new = evaluate_fitness_single(x_new, num_vars, fitness_func_id)
            num_eval = num_eval + 1
            
            IF (f_new >= f + c1 * alpha * directional_deriv) THEN
                RETURN
            END IF
            
            alpha = rho * alpha
        END DO
        
        alpha = 1.0D-6
        DO i = 1, num_vars
            x_new(i) = x(i) + alpha * search_dir(i)
            x_new(i) = DMAX1(lower, DMIN1(upper, x_new(i)))
        END DO
        f_new = evaluate_fitness_single(x_new, num_vars, fitness_func_id)
        num_eval = num_eval + 1
        
    END SUBROUTINE line_search_backtracking
    
    FUNCTION evaluate_fitness_single(x, num_vars, fitness_func_id) RESULT(f)
        USE fitness_functions
        DOUBLE PRECISION, DIMENSION(num_vars), INTENT(IN) :: x
        INTEGER, INTENT(IN) :: num_vars, fitness_func_id
        DOUBLE PRECISION :: f
        
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
            CASE(6)
                f = fitness_slow_rastrigin(x,num_vars)
            CASE DEFAULT
                f = fitness_five_peaks(x, num_vars)
        END SELECT
        
    END FUNCTION evaluate_fitness_single
    
    SUBROUTINE consolidate_converged_niches(centers, center_fits, radius_adaptive, &
                                           num_nich, dims, merge_threshold_factor)
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(INOUT) :: center_fits, radius_adaptive
        INTEGER, INTENT(IN) :: num_nich, dims
        DOUBLE PRECISION, INTENT(IN) :: merge_threshold_factor
        
        INTEGER :: i, j, k, merged_count, survivor_idx, merged_idx
        DOUBLE PRECISION :: dist, fitness_diff, merge_threshold
        LOGICAL :: should_merge
        DOUBLE PRECISION, PARAMETER :: fitness_tol = 0.01D0
        DOUBLE PRECISION, PARAMETER :: invalid_niche = -999.0D0
        
        merged_count = 0
        
        PRINT *, "Consolidating niches that converged to same optimum..."
        PRINT *
        
        DO i = 1, num_nich - 1
            IF (center_fits(i) < -900.0D0) CYCLE
            
            DO j = i + 1, num_nich
                IF (center_fits(j) < -900.0D0) CYCLE
                
                ! Explicit distance calculation
                dist = 0.0D0
                DO k = 1, dims
                    dist = dist + (centers(i,k) - centers(j,k))**2
                END DO
                dist = DSQRT(dist)
                
                fitness_diff = DABS(center_fits(i) - center_fits(j))
                
                merge_threshold = merge_threshold_factor * DMAX1(radius_adaptive(i), radius_adaptive(j))
                
                should_merge = (dist < merge_threshold) .AND. &
                              (fitness_diff < fitness_tol * DMAX1(DABS(center_fits(i)), DABS(center_fits(j))))
                
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
    
    DOUBLE PRECISION, PARAMETER :: invalid_niche = -999.0D0
    DOUBLE PRECISION, PARAMETER :: invalid_check = -998.0D0
    
    ! === CONFIGURABLE PARAMETERS (can be loaded from input file) ===
    ! Variables with default values (will be read from input file later)
    INTEGER :: pop_size = 300
    INTEGER :: max_gen = 50
    INTEGER :: num_vars = 2
    INTEGER :: num_niches_preset = 10 ! preset
    INTEGER :: num_niches  ! Is calculated based on pop_size and num_vars, upper bound
    DOUBLE PRECISION :: min_bound = -40.0D0
    DOUBLE PRECISION :: max_bound = 40.0D0
    DOUBLE PRECISION :: radius_factor = 0.10D0
    DOUBLE PRECISION :: merge_radius_factor = 0.75D0
    INTEGER :: merge_interval = 2
    DOUBLE PRECISION :: consolidate_threshold = 0.2D0
    DOUBLE PRECISION :: cmaes_c_mu = 0.5D0
    DOUBLE PRECISION :: cmaes_damp = 2.0D0
    DOUBLE PRECISION :: diversity_rate = 0.10D0
    INTEGER :: print_interval = 5
    LOGICAL :: print_warnings = .FALSE.
    
    ! These remain as parameters
    INTEGER, PARAMETER :: fitness_function = 1
    DOUBLE PRECISION, PARAMETER :: reg_param = 1.0D-10  ! Tightened for double precision
    
    ! Calculated variables (derived from configurable parameters)
    INTEGER :: min_niche_size  ! Calculated as 2 * num_vars for stable covariance estimation
    
    ! === STATE VARIABLES ===
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cmaes_sigma
    DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: cmaes_cov
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: cmaes_pc, cmaes_ps
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: cmaes_chin
    LOGICAL, DIMENSION(:), ALLOCATABLE :: cmaes_initialized
    
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: population, new_population
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fitness
    INTEGER, DIMENSION(:), ALLOCATABLE :: niche_membership
    LOGICAL, DIMENSION(:), ALLOCATABLE :: niche_converged
    
    DOUBLE PRECISION :: chin_calc
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: niche_centers
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: niche_fitness
    DOUBLE PRECISION, DIMENSION(:, :, :), ALLOCATABLE :: niche_inv_cov
    LOGICAL, DIMENSION(:), ALLOCATABLE :: cov_valid
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: niche_radius_adaptive
    DOUBLE PRECISION :: niche_radius_base
    
    INTEGER :: i,j
    INTEGER :: gen, num_eval, niche_id, member_count
    INTEGER, DIMENSION(:), ALLOCATABLE :: niche_members
    DOUBLE PRECISION :: time_start, time_end
    
    DOUBLE PRECISION, DIMENSION(:, :), ALLOCATABLE :: work_mat
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: work_vec
    
    ! === INITIALIZATION ===
    CALL CPU_TIME(time_start)
    CALL RANDOM_SEED()
    
    ! Calculate derived parameters
    min_niche_size = 2 * num_vars  ! Safe minimum for stable covariance estimation
    
    ! Calculate num_niches to ensure adequate samples per niche for covariance estimation
    num_niches = MIN(num_niches_preset, (pop_size / (3 * num_vars)) - 1)
    
    ! Print configuration
    PRINT *
    PRINT *, "=========================================================="
    PRINT *, "VANESA Configuration"
    PRINT *, "=========================================================="
    PRINT *, "Dimensions:"
    PRINT *, "  Dimensions       = ", num_vars
    PRINT *, "  Population size  = ", pop_size
    PRINT *, "  Retained niches  = ", num_niches, " (calculated: pop_size/(3*num_vars) - 1)"
    PRINT *, "  Maximum gens     = ", max_gen
    PRINT *
    PRINT *, "Search space:"
    PRINT *, "  Min. bound       = ", min_bound
    PRINT *, "  Max. bound       = ", max_bound
    PRINT *
    PRINT *, "Algorithm parameters:"
    PRINT *, "  Min. niche size  = ", min_niche_size, " (= 2 * num_vars)"
    PRINT *, "  Initial radius   = ", radius_factor
    PRINT *, "  Merge radius fac.= ", merge_radius_factor
    PRINT *, "  Merge interval   = ", merge_interval
    PRINT *, "  Diversity inject.= ", diversity_rate
    PRINT *
    PRINT *, "CMA-ES parameters:"
    PRINT *, "  cmaes c_mu       = ", cmaes_c_mu
    PRINT *, "  cmaes damping    = ", cmaes_damp
    PRINT *
    PRINT *, "Printing:"
    PRINT *, "  Print warnings?  =", print_warnings
    PRINT *, "  Print interval   =", print_interval
    PRINT *
    PRINT *, "=========================================================="
    PRINT *
    
    ! Allocate arrays that depend on configurable dimensions
    ALLOCATE(cmaes_sigma(num_niches))
    ALLOCATE(cmaes_cov(num_niches, num_vars, num_vars))
    ALLOCATE(cmaes_pc(num_niches, num_vars))
    ALLOCATE(cmaes_ps(num_niches, num_vars))
    ALLOCATE(cmaes_chin(num_niches))
    ALLOCATE(cmaes_initialized(num_niches))
    ALLOCATE(cov_valid(num_niches))
    ALLOCATE(fitness(pop_size))
    ALLOCATE(new_population(pop_size, num_vars))
    ALLOCATE(niche_membership(pop_size))
    ALLOCATE(niche_centers(num_niches, num_vars))
    ALLOCATE(niche_fitness(num_niches))
    ALLOCATE(niche_inv_cov(num_niches, num_vars, num_vars))
    ALLOCATE(niche_radius_adaptive(num_niches))
    ALLOCATE(niche_members(pop_size))
    ALLOCATE(niche_converged(num_niches))
    ALLOCATE(population(pop_size, num_vars))
    ALLOCATE(work_mat(num_vars, num_vars))
    ALLOCATE(work_vec(num_vars))
    
    ! Scale radius with sqrt(dims) to maintain consistent basin coverage
    ! Normalized to behavior at 2D (test function)
    ! I am assuming the basins have the same volume, which is a reasonable assumption
    niche_radius_base = radius_factor * (max_bound - min_bound) * DSQRT(DBLE(num_vars) / 2.0D0)
    niche_radius_adaptive = niche_radius_base
    niche_converged = .FALSE.
    
    CALL initialize_population_lhs(population, pop_size, num_vars, min_bound, max_bound)
    
    DO i = 1, num_niches
        CALL reset_cmaes_state(i, cmaes_sigma, cmaes_cov, cmaes_pc, cmaes_ps, &
                              cmaes_initialized, num_niches, num_vars)
    END DO
    
    cmaes_chin = DSQRT(chi_square_critical(num_vars, 0.5D0))
    chin_calc = chi_square_critical(num_vars, 0.95D0)
    
    niche_centers = invalid_niche
    niche_fitness = invalid_niche
    cov_valid = .FALSE.
    num_eval = 0
    
    ! === MAIN EVOLUTIONARY LOOP ===
    DO gen = 1, max_gen
        CALL evaluate_fitness(population, fitness, pop_size, num_vars, num_eval)
        
        CALL update_niches_improved(population, fitness, niche_centers, niche_fitness, &
                                    num_niches, pop_size, num_vars, niche_radius_base)
        
        IF (gen == 1) THEN
            CALL assign_niches_euclidean(population, niche_centers, niche_membership, &
                                        pop_size, num_niches, num_vars, niche_radius_adaptive)
        END IF
        
        CALL update_niche_covariances(population, niche_membership, niche_centers, &
                                      niche_inv_cov, cov_valid, &
                                      pop_size, num_niches, num_vars)
        
        CALL assign_niches_mahalanobis(population, niche_centers, niche_inv_cov, &
                                      cov_valid, niche_membership, pop_size, &
                                      num_niches, num_vars, chin_calc)
        
        CALL update_adaptive_niche_radius(population, niche_membership, niche_centers, &
                                         niche_radius_adaptive, niche_radius_base, &
                                         pop_size, num_niches, num_vars)
        
        CALL create_new_generation_cmaes(new_population, niche_membership, &
                                        niche_centers, pop_size, num_vars, num_niches, &
                                        cmaes_sigma, cmaes_cov, min_bound, max_bound)
        
        CALL inject_diversity(new_population, niche_membership, pop_size, num_vars, &
                             min_bound, max_bound, diversity_rate)
        
        population = new_population
        
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
        
        IF (MOD(gen, merge_interval) == 0) THEN
            CALL merge_overlapping_niches(niche_centers, niche_fitness, &
                                         niche_radius_adaptive, cmaes_sigma, &
                                         cmaes_cov, cmaes_pc, cmaes_ps, cmaes_initialized, &
                                         num_niches, num_vars,print_warnings)
            CALL fix_ghost_niches(niche_centers, niche_membership, niche_radius_adaptive, &
                                 cmaes_sigma, pop_size, num_niches, num_vars,print_warnings)
        END IF
        
        IF (MOD(gen, print_interval) == 0) THEN
            CALL print_progress(gen, niche_centers, niche_fitness, niche_radius_adaptive, &
                               cmaes_sigma, num_eval, num_niches, num_vars)
        END IF
    END DO
    
    ! === BFGS REFINEMENT PHASE ===
    PRINT *, "==================================================="
    PRINT *, "Starting BFGS refinement phase"
    PRINT *, "==================================================="
    PRINT *
    
    DO niche_id = 1, num_niches
        IF (is_niche_valid(niche_id, niche_centers, num_niches, num_vars)) THEN
            PRINT *, "Refining niche", niche_id
            
            ! Explicit copy to avoid array slice issues
            DO i = 1, num_vars
                work_vec(i) = niche_centers(niche_id, i)
            END DO
            
            ! Explicit copy of covariance matrix slice
            DO j = 1, num_vars
                DO i = 1, num_vars
                    work_mat(i, j) = cmaes_sigma(niche_id)**2 * cmaes_cov(niche_id, i, j)
                END DO
            END DO
            
            CALL refine_niche_bfgs(work_vec, work_mat, &
                                   niche_fitness(niche_id), num_vars, &
                                   fitness_function, min_bound, max_bound, &
                                   num_eval, niche_converged(niche_id))
            
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
                            cmaes_sigma, niche_converged, num_niches, num_vars)
    
    CALL CPU_TIME(time_end)
    PRINT *, "Time taken: ", time_end - time_start, " seconds"
    PRINT *, "Total evaluations: ", num_eval
    
    ! Deallocate dynamic arrays
    DEALLOCATE(cmaes_sigma)
    DEALLOCATE(cmaes_cov)
    DEALLOCATE(cmaes_pc)
    DEALLOCATE(cmaes_ps)
    DEALLOCATE(cmaes_chin)
    DEALLOCATE(cmaes_initialized)
    DEALLOCATE(cov_valid)
    DEALLOCATE(fitness)
    DEALLOCATE(new_population)
    DEALLOCATE(niche_membership)
    DEALLOCATE(niche_centers)
    DEALLOCATE(niche_fitness)
    DEALLOCATE(niche_inv_cov)
    DEALLOCATE(niche_radius_adaptive)
    DEALLOCATE(niche_members)
    DEALLOCATE(niche_converged)
    DEALLOCATE(population)
    DEALLOCATE(work_mat)
    DEALLOCATE(work_vec)
    
CONTAINS

    ! ==========================================================================
    ! UTILITY FUNCTIONS
    ! ==========================================================================
    
    FUNCTION is_niche_valid(niche_id, centers, num_nich, dims) RESULT(valid)
        INTEGER, INTENT(IN) :: niche_id, num_nich, dims
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        LOGICAL :: valid
        ! Better looking function :-)
        valid = .FALSE.
        IF (niche_id .gt. 0) valid = (centers(niche_id, 1) > invalid_check)

    END FUNCTION is_niche_valid
    
    FUNCTION euclidean_distance_explicit(point1, point2, dims) RESULT(dist)
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: point1, point2
        DOUBLE PRECISION :: dist
        INTEGER :: i
        dist = 0.0D0
        DO i = 1, dims
            dist = dist + (point1(i) - point2(i))**2
        END DO
        dist = DSQRT(dist)
    END FUNCTION euclidean_distance_explicit
    
    FUNCTION count_valid_niches(center_fits, num_nich) RESULT(count)
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(IN) :: center_fits
        INTEGER, INTENT(IN) :: num_nich
        INTEGER :: count, i
        count = 0
        DO i = 1, num_nich
            IF (center_fits(i) > invalid_niche) count = count + 1
        END DO
    END FUNCTION count_valid_niches
    
    SUBROUTINE swap_real(a, b)
        DOUBLE PRECISION, INTENT(INOUT) :: a, b
        DOUBLE PRECISION :: temp
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
        DOUBLE PRECISION, INTENT(IN) :: p
        DOUBLE PRECISION :: z
        DOUBLE PRECISION :: q, r
        DOUBLE PRECISION, PARAMETER :: p_low  = 0.02425D0
        DOUBLE PRECISION, PARAMETER :: p_high = 0.97575D0
        
        DOUBLE PRECISION, PARAMETER :: a(6) = (/ -3.969683028665376D+01,  2.209460984245205D+02, &
                                     -2.759285104469687D+02,  1.383577518672690D+02, &
                                     -3.066479806614716D+01,  2.506628277459239D+00 /)
        DOUBLE PRECISION, PARAMETER :: b(5) = (/ -5.447609879822406D+01,  1.615858368580409D+02, &
                                     -1.556989798598866D+02,  6.680131188771972D+01, &
                                     -1.328068155288572D+01 /)
        
        DOUBLE PRECISION, PARAMETER :: c(6) = (/ -7.784894002430293D-03, -3.223964580411365D-01, &
                                     -2.400758277161838D+00, -2.549732539343734D+00, &
                                      4.374664141464968D+00,  2.938163982698783D+00 /)
        DOUBLE PRECISION, PARAMETER :: d(4) = (/  7.784695709041462D-03,  3.224671290700398D-01, &
                                      2.445134137142996D+00,  3.754408661907416D+00 /)
        
        IF (p < p_low) THEN
            q = DSQRT(-2.0D0 * DLOG(p))
            z = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
                ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0D0)
        ELSE IF (p <= p_high) THEN
            q = p - 0.5D0
            r = q * q
            z = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6)) * q / &
                (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0D0)
        ELSE
            q = DSQRT(-2.0D0 * DLOG(1.0D0 - p))
            z = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
                 ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0D0)
        END IF
        
    END FUNCTION inverse_normal_cdf
    
    FUNCTION chi_square_critical(df, percentile) RESULT(chi_crit)
        INTEGER, INTENT(IN) :: df
        DOUBLE PRECISION, INTENT(IN) :: percentile
        DOUBLE PRECISION :: chi_crit
        DOUBLE PRECISION :: z_alpha
        
        z_alpha = inverse_normal_cdf(percentile)
        chi_crit = df * (1.0D0 - 2.0D0/(9.0D0*df) + z_alpha*DSQRT(2.0D0/(9.0D0*df)))**3
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
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(IN) :: pop
        DOUBLE PRECISION, DIMENSION(size), INTENT(IN) :: fit
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(INOUT) :: center_fits
        INTEGER, INTENT(IN) :: num_nich, size, dims
        DOUBLE PRECISION, INTENT(IN) :: radius
        
        INTEGER :: i, j, k, idx, valid_niches
        INTEGER, DIMENSION(size) :: sorted_indices
        DOUBLE PRECISION :: dist, min_niche_separation
        LOGICAL :: is_new_niche
        DOUBLE PRECISION, DIMENSION(dims) :: temp_point1, temp_point2
        
        CALL sort_by_fitness(fit, sorted_indices, size)
        valid_niches = count_valid_niches(center_fits, num_nich)
        min_niche_separation = radius * 1.00D0
        
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
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(IN) :: pop
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        INTEGER, DIMENSION(size), INTENT(OUT) :: membership
        INTEGER, INTENT(IN) :: size, num_nich, dims
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(IN) :: radius_adaptive
        INTEGER :: i, j, k, closest_niche
        DOUBLE PRECISION :: dist, min_dist
        DOUBLE PRECISION, DIMENSION(dims) :: temp_point1, temp_point2
        
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
            
            IF (closest_niche > 0 .AND. min_dist <= 2.0D0) THEN
                membership(i) = closest_niche
            END IF
        END DO
    END SUBROUTINE assign_niches_euclidean

    SUBROUTINE assign_niches_mahalanobis(pop, centers, inv_cov, cov_valid, membership, &
                                        size, num_nich, dims, chi_threshold)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(IN) :: pop
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(IN) :: inv_cov
        LOGICAL, DIMENSION(num_nich), INTENT(IN) :: cov_valid
        INTEGER, DIMENSION(size), INTENT(OUT) :: membership
        INTEGER, INTENT(IN) :: size, num_nich, dims
        DOUBLE PRECISION, INTENT(IN) :: chi_threshold
        INTEGER :: i, j, k, closest_niche
        DOUBLE PRECISION :: dist, min_dist
        DOUBLE PRECISION, DIMENSION(dims) :: temp_point1, temp_point2
        DOUBLE PRECISION, DIMENSION(dims,dims) :: temp_inv_cov
        
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
                    ! Explicit copy to avoid array slice issues
                    DO k = 1, dims
                        temp_inv_cov(k,1:dims) = inv_cov(j,k,1:dims)
                    END DO
                    dist = mahalanobis_distance_explicit(temp_point1, temp_point2, &
                                                        temp_inv_cov, dims)
                    dist = DSQRT(2.0D0) * dist / DSQRT(chi_threshold)
                ELSE
                    dist = euclidean_distance_explicit(temp_point1, temp_point2, dims)
                END IF
          
                IF (dist < min_dist) THEN
                    min_dist = dist
                    closest_niche = j
                END IF
            END DO
            
            IF (closest_niche > 0 .AND. min_dist <= 2.0D0) THEN
                membership(i) = closest_niche
            END IF
        END DO
    END SUBROUTINE assign_niches_mahalanobis
    
    SUBROUTINE update_adaptive_niche_radius(pop, membership, centers, radius_adaptive, radius_base, &
                                           size, num_nich, dims)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(IN) :: pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(OUT) :: radius_adaptive
        DOUBLE PRECISION, INTENT(IN) :: radius_base
        INTEGER, INTENT(IN) :: size, num_nich, dims
        
        INTEGER :: niche_id, i, j, count
        DOUBLE PRECISION :: mean_dist, dist
        DOUBLE PRECISION, DIMENSION(dims) :: temp_point1, temp_point2
        
        DO niche_id = 1, num_nich
            IF (.NOT. is_niche_valid(niche_id, centers, num_nich, dims)) THEN
                radius_adaptive(niche_id) = radius_base
                CYCLE
            END IF

            count = 0
            mean_dist = 0.0D0
            
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
                mean_dist = mean_dist / DBLE(count)
                radius_adaptive(niche_id) = DMAX1(0.05D0 * radius_base, 1.5D0 * mean_dist)
            END IF
        END DO
    END SUBROUTINE update_adaptive_niche_radius
    
    SUBROUTINE merge_overlapping_niches(centers, center_fits, radius_adaptive, &
                                       sigma, c, pc, ps, initialized, num_nich, dims, print_warn)
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(INOUT) :: center_fits, sigma, radius_adaptive
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: c
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: pc, ps
        LOGICAL, DIMENSION(num_nich), INTENT(INOUT) :: initialized
        INTEGER, INTENT(IN) :: num_nich, dims
        LOGICAL, INTENT(IN) :: print_warn
        
        INTEGER :: i, j, k, merged_count, survivor_idx, merged_idx
        DOUBLE PRECISION :: euclidean_dist, avg_radius, merge_threshold, new_radius
        DOUBLE PRECISION, DIMENSION(dims) :: temp_point1, temp_point2
        
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
                avg_radius = (radius_adaptive(i) + radius_adaptive(j)) / 2.0D0
                merge_threshold = merge_radius_factor * avg_radius
                
                IF (euclidean_dist < merge_threshold) THEN
                    IF (center_fits(i) >= center_fits(j)) THEN
                        survivor_idx = i
                        merged_idx = j
                    ELSE
                        survivor_idx = j
                        merged_idx = i
                    END IF
                    
                    new_radius = DMAX1(radius_adaptive(survivor_idx), &
                                    euclidean_dist + radius_adaptive(merged_idx))
                    radius_adaptive(survivor_idx) = new_radius
                    IF (print_warn) THEN
                        PRINT *, "WARNING: Merged ", i, j, "with centers:", (centers(i,k), k=1,dims), &
                                 " and", (centers(j,k), k=1,dims)
                    ENDIF 
                    CALL invalidate_niche(merged_idx, centers, center_fits, sigma, c, pc, ps, &
                                        initialized, num_nich, dims)
                    
                    merged_count = merged_count + 1
                    IF (survivor_idx /= i) EXIT 
                END IF
            END DO
        END DO
        
        IF ((merged_count > 0).AND.print_warn) THEN
            PRINT *, "WARNING: Merged ", merged_count, " overlapping niches"
        END IF
        
    END SUBROUTINE merge_overlapping_niches
    
    SUBROUTINE fix_ghost_niches(centers, membership, radius_adaptive, sigma, &
                                size, num_nich, dims, print_warn)
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(IN) :: radius_adaptive
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(INOUT) :: sigma
        INTEGER, INTENT(IN) :: size, num_nich, dims
        LOGICAL, INTENT(IN) :: print_warn
        
        INTEGER :: niche_id, member_count, i, k
        DOUBLE PRECISION :: new_sigma
        
        DO niche_id = 1, num_nich
            IF (.NOT. is_niche_valid(niche_id, centers, num_nich, dims)) CYCLE
            
            member_count = 0
            DO i = 1, size
                IF (membership(i) == niche_id) member_count = member_count + 1
            END DO
            
            IF (member_count == 0) THEN
                new_sigma = DMAX1(0.5D0, radius_adaptive(niche_id) * 0.5D0)
                IF (print_warn) THEN
                    PRINT *, "WARNING: Ghost niche ", niche_id, &
                             " at ", (centers(niche_id,k), k=1,dims), &
                             " has 0 members. Resetting sigma from ", sigma(niche_id), &
                             " to ", new_sigma
                ENDIF
                sigma(niche_id) = new_sigma
            END IF
        END DO
    END SUBROUTINE fix_ghost_niches

    SUBROUTINE invalidate_niche(niche_idx, centers, center_fits, sigma, c, pc, ps, &
                               initialized, num_nich, dims)
        INTEGER, INTENT(IN) :: niche_idx, num_nich, dims
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(INOUT) :: center_fits, sigma
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: c
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: pc, ps
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
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(INOUT) :: sigma
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: c
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: pc, ps
        LOGICAL, DIMENSION(num_nich), INTENT(INOUT) :: initialized
        INTEGER :: j, k
        
        sigma(niche_id) = (max_bound - min_bound) / 4.0D0
        DO j = 1, dims
            DO k = 1, dims
                c(niche_id, j, k) = 0.0D0
            END DO
            c(niche_id, j, j) = 1.0D0
            pc(niche_id, j) = 0.0D0
            ps(niche_id, j) = 0.0D0
        END DO
        initialized(niche_id) = .FALSE.
    END SUBROUTINE reset_cmaes_state

    ! ==========================================================================
    ! COVARIANCE AND DISTANCE
    ! ==========================================================================

    FUNCTION mahalanobis_distance_explicit(point1, point2, inv_cov, dims) RESULT(dist)
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: point1, point2
        DOUBLE PRECISION, DIMENSION(dims, dims), INTENT(IN) :: inv_cov
        DOUBLE PRECISION :: dist
        DOUBLE PRECISION, DIMENSION(dims) :: diff, temp
        INTEGER :: i
      
        DO i = 1, dims
            diff(i) = point1(i) - point2(i)
        END DO
        CALL dgemv('N', dims, dims, 1.0D0, inv_cov, dims, diff, 1, 0.0D0, temp, 1)
        dist = DSQRT(DMAX1(0.0D0, DOT_PRODUCT(diff, temp)))
    END FUNCTION mahalanobis_distance_explicit

    SUBROUTINE update_niche_covariances(pop, membership, centers, &
                                      inv_cov_matrices, cov_valid, size, num_nich, dims)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(IN) :: pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(OUT) :: inv_cov_matrices
        LOGICAL, DIMENSION(num_nich), INTENT(OUT) :: cov_valid
        INTEGER, INTENT(IN) :: size, num_nich, dims
        
        INTEGER :: niche_id, i, j, k, member_count, info, idx
        DOUBLE PRECISION, DIMENSION(dims) :: mean_vec, diff
        DOUBLE PRECISION, DIMENSION(dims, dims) :: cov_mat, work_mat
        INTEGER, DIMENSION(size) :: members
        ! LAPACK/BLAS externals
        INTEGER, DIMENSION(dims) :: ipiv
        DOUBLE PRECISION, DIMENSION(4*dims) :: work
        DOUBLE PRECISION :: rcond
        INTEGER, DIMENSION(dims) :: iwork
        
        cov_valid = .FALSE.
        
        DO niche_id = 1, num_nich
            IF (.NOT. is_niche_valid(niche_id, centers, num_nich, dims)) CYCLE
            
            CALL get_niche_members(membership, niche_id, members, member_count, size)
            
            IF (member_count >= 3) THEN
                mean_vec = 0.0D0
                DO i = 1, member_count
                    idx = members(i)
                    DO j = 1, dims
                        mean_vec(j) = mean_vec(j) + pop(idx, j)
                    END DO
                END DO
                DO j = 1, dims
                    mean_vec(j) = mean_vec(j) / DBLE(member_count)
                END DO
                
                cov_mat = 0.0D0
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
                        cov_mat(j, k) = cov_mat(j, k) / DBLE(member_count - 1)
                    END DO
                    cov_mat(j, j) = cov_mat(j, j) + reg_param
                END DO
                
                work_mat = cov_mat
                CALL dgetrf(dims, dims, work_mat, dims, ipiv, info)
                IF (info == 0) THEN
                    CALL dgetri(dims, work_mat, dims, ipiv, work, dims, info)
                    IF (info == 0) THEN
                        CALL dgecon('1', dims, work_mat, dims, 1.0D0, rcond, work, iwork, info)
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
                    CALL dpotrf('U', dims, work_mat, dims, info)
                    IF (info == 0) THEN
                        CALL dpotri('U', dims, work_mat, dims, info)
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
        DOUBLE PRECISION, DIMENSION(size,dims), INTENT(IN) :: population
        DOUBLE PRECISION, DIMENSION(size), INTENT(IN) :: fitness
        INTEGER, DIMENSION(size), INTENT(IN) :: members
        INTEGER, INTENT(IN) :: member_count, dims, niche_id, cov_niche_id
        INTEGER, INTENT(IN) :: pc_niche_id, ps_niche_id, num_nich, size
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: niche_centers
        DOUBLE PRECISION, INTENT(INOUT) :: sigma
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(INOUT) :: cmaes_cov
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(INOUT) :: cmaes_pc, cmaes_ps
        DOUBLE PRECISION, INTENT(IN) :: chin
        LOGICAL, INTENT(INOUT) :: initialized
        
        INTEGER :: i, j, k, idx
        DOUBLE PRECISION, DIMENSION(member_count) :: member_fitness
        DOUBLE PRECISION :: mu_eff, cs, cc, c1, cmu
        DOUBLE PRECISION, DIMENSION(dims) :: mean_step, z_mean, niche_center_copy
        DOUBLE PRECISION, DIMENSION(dims) :: pc_temp, ps_temp
        DOUBLE PRECISION, DIMENSION(dims, dims) :: c_new
        DOUBLE PRECISION :: ps_norm, hsig
    
        
        DO i = 1, member_count
            member_fitness(i) = fitness(members(i))
        END DO
        
        IF (.NOT. initialized) THEN
            initialized = .TRUE.
            RETURN
        END IF
        
        mu_eff = DBLE(MIN(member_count, dims)) / 2.0D0
        cs = (mu_eff + 2.0D0) / (dims + mu_eff + 5.0D0)
        cc = (4.0D0 + mu_eff/dims) / (dims + 4.0D0 + 2.0D0*mu_eff/dims)
        c1 = 2.0D0 / ((dims + 1.3D0)**2 + mu_eff)
        cmu = DMIN1(1.0D0 - c1, 2.0D0*(mu_eff - 2.0D0 + 1.0D0/mu_eff) / ((dims+2.0D0)**2 + mu_eff))
        
        ! Copy niche center
        DO i = 1, dims
            niche_center_copy(i) = niche_centers(niche_id, i)
        END DO
        
        mean_step = 0.0D0
        DO i = 1, MIN(INT(mu_eff), member_count)
            idx = members(i)
            DO j = 1, dims
                mean_step(j) = mean_step(j) + (population(idx, j) - niche_center_copy(j)) / sigma
            END DO
        END DO
        DO j = 1, dims
            mean_step(j) = mean_step(j) / DBLE(MIN(INT(mu_eff), member_count))
        END DO
        
        ! Update ps - Copy to local contiguous array first
        DO i = 1, dims
            ps_temp(i) = cmaes_ps(ps_niche_id, i)
        END DO
        
        DO i = 1, dims
            ps_temp(i) = (1.0D0 - cs) * ps_temp(i) + DSQRT(cs * (2.0D0 - cs)) * mean_step(i)
        END DO
        
        ps_norm = 0.0D0
        DO i = 1, dims
            ps_norm = ps_norm + ps_temp(i)**2
        END DO
        ps_norm = DSQRT(ps_norm)
        
        hsig = 0.0D0
        IF (ps_norm / DSQRT(1.0D0 - (1.0D0-cs)**2) < 1.4D0 * chin) THEN
            hsig = 1.0D0
        END IF
        
        ! Update pc - Copy to local contiguous array first
        DO i = 1, dims
            pc_temp(i) = cmaes_pc(pc_niche_id, i)
        END DO
        
        DO i = 1, dims
            pc_temp(i) = (1.0D0 - cc) * pc_temp(i) + hsig * DSQRT(cc * (2.0D0 - cc)) * mean_step(i)
        END DO
        
        ! Update covariance
        DO k = 1, dims
            DO j = 1, dims
                c_new(j, k) = (1.0D0 - c1 - cmu) * cmaes_cov(cov_niche_id, j, k)
            END DO
        END DO
        
        ! Use local contiguous pc_temp for BLAS call
        CALL dger(dims, dims, c1, pc_temp, 1, pc_temp, 1, c_new, dims)
        
        DO i = 1, MIN(INT(mu_eff), member_count)
            idx = members(i)
            DO j = 1, dims
                z_mean(j) = (population(idx, j) - niche_center_copy(j)) / sigma
            END DO
            CALL dger(dims, dims, cmu/mu_eff, z_mean, 1, z_mean, 1, c_new, dims)
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
        
        sigma = sigma * DEXP((cs / cmaes_damp) * (ps_norm / chin - 1.0D0))
        sigma = DMAX1(1.0D-6 * (max_bound - min_bound), &
                   DMIN1(1.0D0 * (max_bound - min_bound), sigma))
        
    END SUBROUTINE update_cmaes_state_niche
    
    SUBROUTINE create_new_generation_cmaes(new_pop, membership, centers, &
                                           size, dims, num_nich, sigma_cmaes, c_cmaes, &
                                           lower, upper)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(OUT) :: new_pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        INTEGER, INTENT(IN) :: size, dims, num_nich
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(IN) :: sigma_cmaes
        DOUBLE PRECISION, DIMENSION(num_nich, dims, dims), INTENT(IN) :: c_cmaes
        DOUBLE PRECISION, INTENT(IN) :: lower, upper
        
        INTEGER :: i, j, k, niche_id
        DOUBLE PRECISION, DIMENSION(dims) :: sample, center_copy
        DOUBLE PRECISION, DIMENSION(dims, dims) :: cov_copy
        
        DO i = 1, size
            niche_id = membership(i)
            
            IF (niche_id > 0 .AND. is_niche_valid(niche_id, centers, num_nich, dims)) THEN
                ! Explicit copy to ensure contiguous memory
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
                    new_pop(i, j) = DMAX1(lower, DMIN1(upper, sample(j)))
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
        DOUBLE PRECISION, DIMENSION(dims), INTENT(IN) :: mean
        DOUBLE PRECISION, INTENT(IN) :: sigma
        DOUBLE PRECISION, DIMENSION(dims, dims), INTENT(IN) :: c
        DOUBLE PRECISION, DIMENSION(dims), INTENT(OUT) :: sample
        INTEGER, INTENT(IN) :: dims
        
        DOUBLE PRECISION, DIMENSION(dims) :: z, y
        DOUBLE PRECISION, DIMENSION(dims, dims) :: c_chol  ! Cholesky factor
        INTEGER :: i, j, info
        DOUBLE PRECISION :: u1, u2, radius, z1, z2
        DOUBLE PRECISION, PARAMETER :: two_pi = 6.283185307179586D0
        DOUBLE PRECISION, PARAMETER :: max_radius = 5.0D0  ! Clamp extreme values to prevent overflow
        
        ! Generate standard normal random vector z ~ N(0, I)
        DO i = 1, dims, 2
            CALL RANDOM_NUMBER(u1)
            CALL RANDOM_NUMBER(u2)
            
            ! Protect against extreme values that can cause overflow
            u1 = DMAX1(u1, 1.0D-8)  ! Ensure u1 is not too close to 0
            
            ! Box-Muller transform with clamping
            radius = DSQRT(-2.0D0 * DLOG(u1))
            radius = DMIN1(radius, max_radius)  ! Prevent extreme outliers
            
            z1 = radius * DCOS(two_pi * u2)
            z2 = radius * DSIN(two_pi * u2)
            
            z(i) = z1
            IF (i + 1 <= dims) THEN
                z(i + 1) = z2
            END IF
        END DO
        
        ! Copy covariance matrix for Cholesky decomposition
        DO j = 1, dims
            DO i = 1, dims
                c_chol(i, j) = c(i, j)
            END DO
        END DO
        
        ! Perform Cholesky decomposition c = L * L^T
        ! This is the mathematically correct way to sample from multivariate normal
        ! We need to use L (the Cholesky factor), not C (the full covariance)
        CALL dpotrf('L', dims, c_chol, dims, info)
        
        IF (info /= 0) THEN
            ! Cholesky failed - matrix not positive definite
            ! Fall back to identity matrix (isotropic sampling)
            PRINT *, "WARNING: Cholesky decomposition failed (info=", info, &
                     "), using isotropic sampling"
            DO i = 1, dims
                DO j = 1, dims
                    IF (i == j) THEN
                        c_chol(i, j) = 1.0D0
                    ELSE
                        c_chol(i, j) = 0.0D0
                    END IF
                END DO
            END DO
        ELSE
            ! Zero out upper triangle (dpotrf only fills lower triangle)
            DO i = 1, dims
                DO j = i+1, dims
                    c_chol(i, j) = 0.0D0
                END DO
            END DO
        END IF
        
        ! Compute y = L * z (using Cholesky factor L, not full covariance C)
        ! This is correct: x = mean + sigma * L * z where C = L * L^T
        CALL dgemv('N', dims, dims, 1.0D0, c_chol, dims, z, 1, 0.0D0, y, 1)
        
        ! Final sample
        DO i = 1, dims
            sample(i) = mean(i) + sigma * y(i)
        END DO
    
    END SUBROUTINE sample_from_cmaes_distribution  

    ! ==========================================================================
    ! POPULATION MANAGEMENT
    ! ==========================================================================
    
    SUBROUTINE initialize_population_lhs(pop, size, dims, lower, upper)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(OUT) :: pop
        INTEGER, INTENT(IN) :: size, dims
        DOUBLE PRECISION, INTENT(IN) :: lower, upper
        INTEGER :: i, j, k
        DOUBLE PRECISION :: r
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
                pop(i, j) = lower + (upper - lower) * (perm(i) - 1 + r) / DBLE(size)
            END DO
        END DO
    END SUBROUTINE initialize_population_lhs
    
    SUBROUTINE inject_diversity(pop, membership, size, dims, lower, upper, rate)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(INOUT) :: pop
        INTEGER, DIMENSION(size), INTENT(IN) :: membership
        INTEGER, INTENT(IN) :: size, dims
        DOUBLE PRECISION, INTENT(IN) :: lower, upper, rate
        
        INTEGER :: i, j, num_inject
        DOUBLE PRECISION :: r
        DOUBLE PRECISION, DIMENSION(dims) :: temp_vec
        
        num_inject = INT(rate * DBLE(size))
        
        DO i = 1, size
            IF (i <= num_inject) THEN
                CALL random_uniform_vector(temp_vec, dims, lower, upper)
                DO j = 1, dims
                    pop(i, j) = temp_vec(j)
                END DO
            ELSE IF (membership(i) == 0) THEN
                CALL RANDOM_NUMBER(r)
                IF (r < 0.3D0) THEN
                    CALL random_uniform_vector(temp_vec, dims, lower, upper)
                    DO j = 1, dims
                        pop(i, j) = temp_vec(j)
                    END DO
                END IF
            END IF
        END DO
    END SUBROUTINE inject_diversity
    
    SUBROUTINE random_uniform_vector(vec, dims, lower, upper)
        DOUBLE PRECISION, DIMENSION(dims), INTENT(OUT) :: vec
        INTEGER, INTENT(IN) :: dims
        DOUBLE PRECISION, INTENT(IN) :: lower, upper
        INTEGER :: i
        DOUBLE PRECISION :: r
        
        DO i = 1, dims
            CALL RANDOM_NUMBER(r)
            vec(i) = lower + r * (upper - lower)
        END DO
    END SUBROUTINE random_uniform_vector
    
    ! ==========================================================================
    ! FITNESS EVALUATION AND SORTING
    ! ==========================================================================
    
    SUBROUTINE evaluate_fitness(pop, fit, size, dims, evalfs)
        DOUBLE PRECISION, DIMENSION(size, dims), INTENT(IN) :: pop
        DOUBLE PRECISION, DIMENSION(size), INTENT(OUT) :: fit
        INTEGER, INTENT(IN) :: size, dims
        INTEGER, INTENT(INOUT) :: evalfs
        INTEGER :: i, j
        DOUBLE PRECISION, DIMENSION(dims) :: temp_point
        
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, temp_point) SCHEDULE(DYNAMIC)
        DO i = 1, size
            
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
                CASE(6)
                    fit(i) = fitness_slow_rastrigin(temp_point, dims)
                CASE DEFAULT
                    fit(i) = fitness_five_peaks(temp_point, dims)
            END SELECT
        END DO
        !$OMP END PARALLEL DO
        
        ! Update evaluation counter after parallel region
        evalfs = evalfs + size
    END SUBROUTINE evaluate_fitness
    
    SUBROUTINE sort_by_fitness(fitness, indices, size)
        INTEGER, INTENT(IN) :: size
        DOUBLE PRECISION, DIMENSION(size), INTENT(IN) :: fitness
        INTEGER, DIMENSION(size), INTENT(OUT) :: indices
        
        DOUBLE PRECISION, DIMENSION(size) :: fitness_copy
        INTEGER :: i
        
        DO i = 1, size
            indices(i) = i
            fitness_copy(i) = fitness(i)
        END DO
        
        CALL quicksort_descending(fitness_copy, indices, 1, size)
        
    END SUBROUTINE sort_by_fitness

    RECURSIVE SUBROUTINE quicksort_descending(fitness, indices, left, right)
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: fitness
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
        DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: fitness
        INTEGER, DIMENSION(:), INTENT(INOUT) :: indices
        INTEGER, INTENT(IN) :: left, right
        INTEGER, INTENT(OUT) :: pivot_idx
        
        DOUBLE PRECISION :: pivot_value
        INTEGER :: i
        
        pivot_idx = (left + right) / 2
        pivot_value = fitness(pivot_idx)
        
        CALL swap_real(fitness(pivot_idx), fitness(right))
        CALL swap_int(indices(pivot_idx), indices(right))
        
        pivot_idx = left
        DO i = left, right - 1
            IF (fitness(i) > pivot_value .OR. &
                (DABS(fitness(i) - pivot_value) < 1.0D-10 .AND. indices(i) < indices(right))) THEN
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
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(IN) :: center_fits, radius_adaptive, sigma
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
    
    SUBROUTINE print_final_results(centers, center_fits, radius_adaptive, sigma, converged, num_nich, dims)
        DOUBLE PRECISION, DIMENSION(num_nich, dims), INTENT(IN) :: centers
        DOUBLE PRECISION, DIMENSION(num_nich), INTENT(IN) :: center_fits, radius_adaptive, sigma
        INTEGER, INTENT(IN) :: num_nich, dims
        LOGICAL, DIMENSION(num_nich), INTENT(IN) :: converged
        
        INTEGER :: i, j, k, valid_count, unconverged_count
        INTEGER, DIMENSION(num_nich) :: sorted_indices
        DOUBLE PRECISION, DIMENSION(num_nich) :: temp_fitness
        
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
                IF (.NOT. converged(i)) unconverged_count = unconverged_count + 1
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
            IF (.NOT. converged(i)) PRINT *, "** WARNING: Unconverged niche (BFGS) - check the log! **" 
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
