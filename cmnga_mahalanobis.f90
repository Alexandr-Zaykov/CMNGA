program CumulativeMultiNichingGA
    implicit none
    
    ! External BLAS/LAPACK functions
    external :: sgemv, sger, sgetrf, sgetri, sgecon, spotrf, spotri
    ! chi_square tables
    ! For Mahalanobis distance, use chi-square threshold
    ! Higher than 30D is very likely completely useless, but the code snippet may help someone
    ! Chi-square threshold for Mahalanobis distance
    ! It provides a multivariate test whether a point is an outlier = hypothesis
    ! If we are following statistics, 95% threshold should be used
    ! However, this provides VERY broad niches... unless sampled enough...
    ! We will likely deal with rather poor sampling and might need to use broader definition of an outlier.
    ! = large p-value
    ! NOTE: This still behaves a little counter-intuitively...
    !       It may also not be the best idea and we might want to go with some empirical constant.
    real,parameter :: chi_table_01(37)=(/10.828, 13.816, 16.266, 18.467, 20.515, &
                                 & 22.458, 24.322, 26.124, 27.877, 29.588, &
                                 & 31.264, 32.909, 34.528, 36.123, 37.697, &
                                 & 39.252, 40.790, 42.312, 43.820, 45.315, &
                                 & 46.797, 48.268, 49.728, 51.179, 52.620, &
                                 & 54.052, 55.476, 56.892, 58.302, 59.703, &
                                 & 73.402, 86.661, 99.607,112.317,124.839, &
                                 & 137.208,149.449/) ! At 0.1% significance
    real,parameter :: chi_table_1(37)=(/ 6.635,  9.210, 11.345, 13.277, 15.086, &
                                 &16.812, 18.475, 20.090, 21.666, 23.209, &
                                 &24.725, 26.217, 27.688, 29.141, 30.578, &
                                 &32.000, 33.409, 34.805, 36.191, 37.566, &
                                 &38.932, 40.289, 41.638, 42.980, 44.314, &
                                 &45.642, 46.963, 48.278, 49.588, 50.892, &
                                 &63.691, 76.154, 88.379,100.425,112.329, &
                                 &124.116,135.807/) ! At 1% significance
    real,parameter :: chi_table_5(37)=(/ 3.841, 5.991, 7.815, 9.488,11.070,12.592,14.067,15.507,16.919,18.307,&
                                     &19.675,21.026,22.362,23.685,24.996,26.296,27.587,28.869,30.144,31.410,&
                                     &32.671,33.924,35.172,36.415,37.652,38.885,40.113,41.337,42.557,43.773,&
                                     &55.758,67.505,79.082,90.531,101.879,113.145,124.342/) ! At 5% significance
    real,parameter :: chi_table_10(37)=(/ 2.706, 4.605, 6.251, 7.779, 9.236,10.645,12.017,13.362,14.684,15.987,&
                                        &17.275,18.549,19.812,21.064,22.307,23.542,24.769,25.989,27.204,28.412,&
                                        &29.615,30.813,32.007,33.196,34.382,35.563,36.741,37.916,39.087,40.256,&
                                        &51.805,63.167,74.397,85.527, 96.578,107.565,118.498/) ! At 10% significance
    ! Bad idea how to provide tighter niches:
    real,parameter :: chi_table_50(37)=(/ 0.455, 1.386, 2.366, 3.357, 4.351, &
                                & 5.348, 6.346, 7.344, 8.343, 9.342, &
                                &10.341,11.340,12.340,13.339,14.339, &
                                &15.338,16.338,17.338,18.338,19.337, &
                                &20.337,21.337,22.337,23.337,24.336, &
                                &25.336,26.336,27.336,28.336,29.336, &
                                &39.335,49.335,59.335,69.334,79.334, &
                                &89.334,99.334/) ! At 50% significance
    real,parameter :: chi_table_70(37)=(/ 0.148, 0.713, 1.424, 2.195, 2.999, &
                                & 3.828, 4.671, 5.527, 6.393, 7.267, &
                                & 8.148, 9.034, 9.926,10.821,11.721, &
                                &12.624,13.531,14.440,15.352,16.266, &
                                &17.182,18.101,19.021,19.943,20.867, &
                                &21.792,22.719,23.647,24.577,25.508, &
                                &35.534,45.616,55.758,65.950,76.188, &
                                &86.464,96.774/) ! At 70% significance
    real,parameter :: chi_table_90(37)=(/ 0.016, 0.211, 0.584, 1.064, 1.610, &
                                & 2.204, 2.833, 3.490, 4.168, 4.865, &
                                & 5.578, 6.304, 7.042, 7.790, 8.547, &
                                & 9.312,10.085,10.865,11.651,12.443, &
                                &13.240,14.041,14.848,15.659,16.473, &
                                &17.292,18.114,18.939,19.768,20.599, &
                                &29.588,39.364,49.475,59.770,70.195, &
                                &80.725,91.342/) ! At 90% significance
    real,parameter :: chi_table_95(37)=(/ 0.004, 0.103, 0.352, 0.711, 1.145, &
                                & 1.635, 2.167, 2.733, 3.325, 3.940, &
                                & 4.575, 5.226, 5.892, 6.571, 7.261, &
                                & 7.962, 8.672, 9.390,10.117,10.851, &
                                &11.591,12.338,13.091,13.848,14.611, &
                                &15.379,16.151,16.928,17.708,18.493, &
                                &26.509,35.249,44.314,53.540,62.830, &
                                &72.153,81.500/)
    
    ! Parameters
    ! There are further distribution parameters for crossover and mutation "eta_c" – we use 20 for both for now.
    integer, parameter :: pop_size = 50 ! Population size
    integer, parameter :: max_gen = 400    ! Maximum generations
    integer, parameter :: num_vars = 1           ! Number of variables in optimization (up to 30D)
    integer, parameter :: num_niches = 20        ! Number of niches to maintain
    integer, parameter :: min_niche_size = 3         ! Minimum current individuals per niche for covariance
    real, parameter :: mutate_rate = 0.25         ! Mutation rate
    real, parameter :: cross_rate = 0.65           ! Crossover rate
    real, parameter :: radius_factor = 0.05       ! Niche radius factor
    real, parameter :: min_bound = -4.0           ! Lower bound for variables
    real, parameter :: max_bound = 4.0            ! Upper bound for variables
    real, parameter :: reg_param = 1e-6           ! Regularization parameter for covariance
  ! New parameters for adaptive memory
    real, parameter :: alpha_start = 0.95   ! High weight for current gen early on (exploration)
    real, parameter :: alpha_end = 0.3      ! Lower weight later (exploitation/stability)
    real, parameter :: niche_stability_threshold = 0.01   ! Distance threshold for "stable" niche
    integer, parameter :: transition_gens = 90*max_gen/100 + 1  ! Generations over which to transition – arbitrary. +1 to avoid acc. 0.
  ! New arrays to track niche history
    real, dimension(num_niches, num_vars) :: prev_niche_centers
    real, dimension(num_niches) :: niche_alpha  ! Individual alpha per niche
    integer, dimension(num_niches) :: niche_stability_count
    logical, dimension(num_niches) :: cov_initialized=.false., niche_stable=.false.
    
    ! Variables
    real, dimension(pop_size, num_vars) :: population      ! Current population
    real, dimension(pop_size) :: fitness                   ! Fitness values
    real, dimension(pop_size, num_vars) :: new_population  ! Next generation
    real, dimension(num_niches, num_vars) :: niche_centers ! Centers of niches
    real, dimension(num_niches) :: niche_fitness           ! Fitness of niche centers
    integer, dimension(pop_size) :: niche_membership       ! Which niche each individual belongs to
    real, dimension(num_niches, num_vars, num_vars) :: niche_covariance ! Covariance matrices
    real, dimension(num_niches, num_vars, num_vars) :: niche_inv_cov    ! Inverse covariance matrices
    logical, dimension(num_niches) :: cov_valid            ! Whether covariance is valid
    real :: niche_radius                                   ! Radius for niching
    integer :: i, gen, num_eval
    real :: time_start, time_end
    
    ! Timing
    call CPU_TIME(time_start)
    ! Initialize random seed
    call random_seed()
    
    ! Set niche radius as a fraction of the search space
    niche_radius = radius_factor * (max_bound - min_bound)
    
    ! Initialize population randomly
    call initialize_population(population, pop_size, num_vars, min_bound, max_bound)
    
    ! Initialize niche centers to invalid values
    niche_centers = -999.0
    niche_fitness = -999.0
    cov_valid = .false.
    ! Debug
    num_eval = 0
    
    ! Main loop
    do gen = 1, max_gen
        ! Evaluate fitness of current population
        call evaluate_fitness(population, fitness, pop_size, num_vars, num_eval)
        
        ! Update niche centers
       !call update_niches(population, fitness, niche_centers, niche_fitness, &
       !                  &num_niches, pop_size, num_vars, niche_radius)
        call update_niches_improved(population, fitness, niche_centers, niche_fitness, &
                           num_niches, pop_size, num_vars, niche_radius)

      ! This just does not work at all.
      ! call update_niches_mahalanobis(population, fitness, niche_centers, niche_fitness, &
      !                               &num_niches, pop_size, num_vars, niche_radius,cov_valid,&
      !                                &niche_inv_cov,chi_table_01(num_vars))
        
        ! Assign individuals to niches
        ! Legacy code. Curiously, it adds some niches. Provides... better results?
        call assign_niches(population, niche_centers, niche_membership, &
                           &pop_size, num_niches, num_vars, niche_radius)
        ! Assign individuals using Mahalanobis distance where possible
        ! Falls back to Euclid whereever needed
        ! Provides... worse results?
        ! Often has problems with few niches.
      ! call assign_niches_mahalanobis(population, niche_centers, niche_inv_cov, &
      !                                &cov_valid, niche_membership, pop_size, &
      !                                &num_niches, num_vars, niche_radius,chi_table_5(num_vars))
        
        call update_adaptive_alpha(niche_centers, prev_niche_centers, niche_stability_count, &
                                  &niche_alpha, gen, num_niches, num_vars, niche_stability_threshold,&
                                  &alpha_start, alpha_end, niche_stable)
        ! Update covariance matrices for each niche
        call update_niche_covariances(population, niche_membership, niche_centers, &
                                     &niche_covariance, niche_inv_cov, cov_valid, &
                                     &pop_size, num_niches, num_vars,&
                                     &niche_alpha,cov_initialized)
        
        ! Re-assign individuals using Mahalanobis distance where possible
        call assign_niches_mahalanobis(population, niche_centers, niche_inv_cov, &
                                      &cov_valid, niche_membership, pop_size, &
                                      &num_niches, num_vars, niche_radius,chi_table_5(num_vars))
      ! ! Update niche centers
        
        ! Create new population
        call create_new_generation(population, fitness, new_population, niche_membership, &
                                  &pop_size, num_vars, num_niches, cross_rate, mutate_rate, &
                                  &min_bound, max_bound)
        
        ! Replace old population with new one
        population = new_population
        
        ! Print progress
        if (mod(gen, 10) == 0) then
            print *, "Generation ", gen
            print *, "Best niches found:"
            do i = 1, num_niches
                if (niche_fitness(i) > -999.0) then
                    if (cov_valid(i)) then 
                      print *, "Niche ", i, ": Position = ", niche_centers(i,:), &
                           " Fitness = ", niche_fitness(i), " Mahalanobis;", " Stable?", niche_stable(i)
                    else
                      print *, "Niche ", i, ": Position = ", niche_centers(i,:), &
                           " Fitness = ", niche_fitness(i), " Euclidean;","   Stable?", niche_stable(i)
                    endif
                end if
            end do
            print *
        end if
    end do
    
    ! Print final results
    print *, "Optimization completed"
    print *, "Final niches found:"
    do i = 1, num_niches
        if (niche_fitness(i) > -999.0) then
          if (cov_valid(i)) then 
            print *, "Niche ", i, ": Position = ", niche_centers(i,:), &
                 " Fitness = ", niche_fitness(i)-1, " Mahalanobis;", " Stable?", niche_stable(i)
          else
            print *, "Niche ", i, ": Position = ", niche_centers(i,:), &
                 " Fitness = ", niche_fitness(i)-1, " Euclidean;", "   Stable?", niche_stable(i)
          endif
        end if
    end do
    call CPU_TIME(time_end)
    print *, "Time taken: ", time_end - time_start, " s"
    print *, "Anticipated Results: +-0.31,0.93,1.55,2.17,2.79,3.41"
    print *, "Number of evaluations: ", num_eval
contains

    function det_posdef(A, n) result(det)
    implicit none
    
    ! Arguments
    integer, intent(in) :: n
    real, intent(in) :: A(n,n)
    real :: det
    
    ! Local variables
    real :: A_copy(n,n)
    integer :: info, i
    
    ! External LAPACK routine
    external :: spotrf
    
    ! Copy the matrix since DPOTRF modifies it
    A_copy = A
    
    ! Compute Cholesky factorization: A = U^T * U
    call spotrf('U', n, A_copy, n, info)
    
    ! Check if factorization was successful
    if (info /= 0) then
        if (info < 0) then
            write(*,*) 'Error: Illegal argument in DPOTRF at position', -info
            call exit(8)
        else
            write(*,*) 'Error: Matrix is not positive definite'
            call exit(8)
        end if
        det = 0.0d0
        return
    end if
    
    ! Compute determinant as product of squares of diagonal elements
    ! det(A) = det(U^T) * det(U) = [det(U)]^2 = [prod(U_ii)]^2
    det = 1.0d0
    do i = 1, n
        det = det * A_copy(i,i)**2
    end do
    
    end function det_posdef


  ! Alternative version using LU decomposition (more general but less efficient for pos def matrices)
    function det_general(A, n) result(det)
    implicit none
    
    ! Arguments
    integer, intent(in) :: n
    real, intent(in) :: A(n,n)
    real :: det
    
    ! Local variables
    real :: A_copy(n,n)
    integer :: ipiv(n), info, i
    
    ! External LAPACK routine
    external :: sgetrf
    
    ! Copy the matrix
    A_copy = A
    
    ! Compute LU factorization with partial pivoting
    call sgetrf(n, n, A_copy, n, ipiv, info)
    
    if (info /= 0) then
        write(*,*) 'Error in LU factorization'
        det = 0.0d0
        return
    end if
    
    ! Compute determinant from diagonal elements and pivot count
    det = 1.0d0
    do i = 1, n
        if (ipiv(i) /= i) then
            det = -det  ! Account for row swaps
        end if
        det = det * A_copy(i,i)
    end do
    
    end function det_general
    ! Initialize population randomly
    subroutine initialize_population(pop, size, dims, lower, upper)
        real, dimension(size, dims), intent(out) :: pop
        integer, intent(in) :: size, dims
        real, intent(in) :: lower, upper
        integer :: i, j
        real :: r
        
        do i = 1, size
            do j = 1, dims
                call random_number(r)
                pop(i, j) = lower + r * (upper - lower)
            end do
        end do
    end subroutine initialize_population
    
    ! Evaluate fitness (using a multi-modal test function)
    subroutine evaluate_fitness(pop, fit, size, dims,evalfs)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(size), intent(out) :: fit
        integer, intent(in) :: size, dims
        integer, intent(inout) :: evalfs
        integer :: i,j

        real                               :: A(5), B(5), H(5), W(5), F = 0.0
 
        data A /-20.0, -5.0, 0.0, 30.0, 30.0/
        data B /-20.0, -25.0, 30.0, 0.0, -30.0/
        data H /0.4, 0.2, 0.7, 1.0, 0.05/
        data W /0.02, 0.5, 0.01, 2.0, 0.1/

       
        ! ADD THE SEARCH FUNCTION HERE
        do i = 1, size
            evalfs = evalfs + 1
            ! Example multi-modal function: Sum of sines
            if (dims == 1) then
                  fit(i) = sin(5.0 * pop(i,1))**6 * exp(-pop(i,1)**2)
                
                ! Results: +-0.31,0.93,1.55,2.17,2.79,3.41
            else
                ! 2D or higher: Combination of sines (only use first 2 dimensions for test)
                ! fit(i) = sin(5.0 * pop(i,1))**2 * sin(5.0 * pop(i,2))**2 * &
                !        exp(-(pop(i,1)**2 + pop(i,2)**2)/2)
                ! Different test function
                F = 0.0 
                do j = 1,5
                    F = F + H(j) / (1.0 + W(j)*( (pop(i,1)-A(j))**2 + (pop(i,2)-B(j))**2 ))
                enddo
                fit(i) = F

                ! Add contribution from other dimensions if present
                if (dims > 2) then
                    do j = 3, dims
                        fit(i) = fit(i) * exp(-(pop(i,j)**2)/4)  ! Decay for higher dims
                    end do
                end if
            end if
            fit(i) = fit(i) + 1.0
        end do
    end subroutine evaluate_fitness
    
    ! Calculate Mahalanobis distance
    function mahalanobis_distance(point1, point2, inv_cov, dims) result(dist)
        integer, intent(in) :: dims
        real, dimension(dims), intent(in) :: point1, point2
        real, dimension(dims, dims), intent(in) :: inv_cov
        real :: dist
        real, dimension(dims) :: diff
        real, dimension(dims) :: temp
       
        diff = point1 - point2
        ! temp = inv_cov * diff (using BLAS SGEMV for efficiency)
        call sgemv('N', dims, dims, 1.0, inv_cov, dims, diff, 1, 0.0, temp, 1)
        ! dist = sqrt(diff^T * temp) (using intrinsic – I failed to use sdot)
        dist = sqrt(max(0.0, dot_product(diff,temp)))
    end function mahalanobis_distance

    subroutine update_adaptive_alpha(centers, prev_centers, stability_count, alpha, &
                                generation, num_nich, dims, stab_threshold,alpha_start,alpha_end,stable)
    real, dimension(num_nich, dims), intent(in) :: centers
    real, dimension(num_nich, dims), intent(inout) :: prev_centers
    integer, dimension(num_nich), intent(inout) :: stability_count
    real, dimension(num_nich), intent(out) :: alpha
    logical,dimension(num_nich), intent(out) :: stable
    integer, intent(in) :: generation, num_nich, dims
    real, intent(in) :: stab_threshold,alpha_start,alpha_end
    
    integer :: i, j
    real :: base_alpha, stability_factor, distance, max_stable_gens
    real, parameter :: stability_bonus = 0.2  ! How much to reduce alpha for stable niches
    
    ! Calculate base alpha that changes with generation
    if (generation <= transition_gens) then
        ! Linear transition from alpha_start to alpha_end
        base_alpha = alpha_start - (alpha_start - alpha_end) * &
                    real(generation - 1) / real(transition_gens)
    else
        base_alpha = alpha_end
    end if
    
    ! Adjust alpha for each niche based on stability
    do i = 1, num_nich
        if (centers(i, 1) > -998.0) then  ! Valid niche
            
            ! Calculate how much this niche has moved
            distance = 0.0
            if (prev_centers(i, 1) > -998.0) then  ! Had previous position
                do j = 1, dims
                    distance = distance + (centers(i, j) - prev_centers(i, j))**2
                end do
                distance = sqrt(distance)
                
                ! Update stability count
                if (distance < stab_threshold) then
                    stability_count(i) = stability_count(i) + 1
                else
                    stability_count(i) = 0  ! Reset if niche moved significantly
                end if
            else
                stability_count(i) = 0  ! New niche
            end if
            
            ! Calculate stability factor (0 = unstable, 1 = very stable)
            max_stable_gens = min(50.0, real(generation))  ! Max stability window, never 0. At least 50 gens.
            stability_factor = min(1.0, real(stability_count(i)) / max_stable_gens)
            
            ! Reduce alpha for stable niches (more memory)
            ! we use the factor^6 for slow ramp up toward 1 in the beginning and focusing in the end
            alpha(i) = base_alpha - stability_bonus * (stability_factor**6)
            ! Check for the stability of the niche
            IF(stability_factor.GT.0.99) stable(i)=.true.
            alpha(i) = max(0.1, min(0.95, alpha(i)))  ! Clamp to reasonable range
            
            ! Update previous center for next iteration
            prev_centers(i, :) = centers(i, :)
        else
            niche_alpha(i) = base_alpha  ! Default for invalid niches
        end if
    end do
    end subroutine update_adaptive_alpha

    ! Update covariance matrices for each niche using LAPACK
    ! ADAPTIVE
    subroutine update_niche_covariances(pop, membership, centers, cov_matrices, &
                                      &inv_cov_matrices, cov_valid, size, num_nich, dims,&
                                      &alpha,cov_init)
        real, dimension(size, dims), intent(in) :: pop
        integer, dimension(size), intent(in) :: membership
        real, dimension(num_nich, dims), intent(in) :: centers
        real, dimension(num_nich, dims, dims), intent(out) :: cov_matrices
        real, dimension(num_nich, dims, dims), intent(out) :: inv_cov_matrices
        logical, dimension(num_nich), intent(out) :: cov_valid
        integer, intent(in) :: size, num_nich, dims
        real, dimension(num_nich), intent(in) :: alpha
        logical, dimension(num_nich), intent(inout) :: cov_init
        
        integer :: niche_id, i, j, k, count, info
        real, dimension(dims) :: mean_vec, diff
        real, dimension(dims, dims) :: cov_mat, work_mat, blended_cov
        integer, dimension(size) :: niche_members
        integer, dimension(dims) :: ipiv
        real, dimension(dims) :: work
        real :: rcond
        integer, dimension(dims) :: iwork
        
        ! External LAPACK functions
        external :: sgetrf, sgetri, sgecon
        
        cov_valid = .false.
        
        do niche_id = 1, num_nich
            if (centers(niche_id, 1) < -998.0) cycle  ! Skip invalid niches
            
            ! Count members in this niche
            count = 0
            do i = 1, size
                if (membership(i) == niche_id) then
                    count = count + 1
                    niche_members(count) = i
                end if
            end do
            
            ! Need at least min_niche_size individuals for reliable covariance
            if (count >= min_niche_size) then
                ! Calculate mean
                mean_vec = 0.0
                do i = 1, count
                    mean_vec = mean_vec + pop(niche_members(i), :)
                end do
                mean_vec = mean_vec / real(count)
                
                ! Calculate covariance matrix using BLAS for efficiency
                cov_mat = 0.0
                do i = 1, count
                    diff = pop(niche_members(i), :) - mean_vec
                    ! Outer product: cov_mat += diff * diff^T
                    call sger(dims, dims, 1.0, diff, 1, diff, 1, cov_mat, dims)
                end do
                cov_mat = cov_mat / real(count - 1)
                
                ! Add regularization to diagonal to ensure positive definiteness
                do j = 1, dims
                    cov_mat(j, j) = cov_mat(j, j) + reg_param
                end do
                
                  ! Adaptive blending with previous covariance
                if (.not. cov_init(niche_id)) then
                    ! First time - use current calculation
                    blended_cov = cov_mat
                    cov_init(niche_id) = .true.
                else
                    ! Blend with previous covariance using adaptive alpha
                    blended_cov = alpha(niche_id) * cov_mat + &
                                 (1.0 - alpha(niche_id)) * cov_matrices(niche_id, :, :)
                end if
                ! Store covariance matrix
                cov_matrices(niche_id, :, :) = blended_cov
                
                ! Calculate inverse using LAPACK
                work_mat = cov_mat  ! Copy since LAPACK destroys input
                ! LU factorization
                call sgetrf(dims, dims, work_mat, dims, ipiv, info)
                
                if (info == 0) then
                    ! Check condition number to avoid nearly singular matrices
                    call sgecon('1', dims, work_mat, dims, 1.0, rcond, work, iwork, info)
                    
                    if (info == 0 .and. rcond > 1e-12) then
                        ! Matrix is well-conditioned, compute inverse
                        call sgetri(dims, work_mat, dims, ipiv, work, dims, info)
                        
                        if (info == 0) then
                            inv_cov_matrices(niche_id, :, :) = work_mat
                            cov_valid(niche_id) = .true.
                        end if
                    end if
                end if
                
                ! If LAPACK inversion failed, try Cholesky decomposition for SPD matrices
                ! If Cholesky fails, default to Euclidean norm 
                if (.not. cov_valid(niche_id)) then
                    work_mat = cov_mat
                    call spotrf('U', dims, work_mat, dims, info)
                    if (info == 0) then
                        call spotri('U', dims, work_mat, dims, info)
                        if (info == 0) then
                            ! Copy upper triangle to lower triangle
                            do j = 1, dims
                                do k = j+1, dims
                                    work_mat(k, j) = work_mat(j, k)
                                end do
                            end do
                            inv_cov_matrices(niche_id, :, :) = work_mat
                            cov_valid(niche_id) = .true.
                        end if
                    end if
                end if
            end if
        end do
    end subroutine update_niche_covariances

        ! Update niche centers based on current population
    subroutine update_niches(pop, fit, centers, center_fits, num_nich, size, dims, radius)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(size), intent(in) :: fit
        real, dimension(num_nich, dims), intent(inout) :: centers
        real, dimension(num_nich), intent(inout) :: center_fits
        integer, intent(in) :: num_nich, size, dims
        real, intent(in) :: radius
        integer :: i, j, k, best_idx
        real :: best_fit, dist
        logical :: is_new_niche
        
        ! First, find the best individual
        best_fit = -1.0
        best_idx = 1
        do i = 1, size
            if (fit(i) > best_fit) then
                best_fit = fit(i)
                best_idx = i
            end if
        end do
        
        ! Check if we need to add a new niche center
        if (center_fits(1) < 0.0) then
            ! First niche - just add it
            centers(1,:) = pop(best_idx,:)
            center_fits(1) = fit(best_idx)
        else
            ! Check if the best individual is close to any existing niche
            is_new_niche = .true.
            do j = 1, num_nich
                if (center_fits(j) < 0.0) exit  ! No more valid niches
                
                ! Calculate Euclidean distance to this niche center
                dist = 0.0
                do k = 1, dims
                    dist = dist + (pop(best_idx,k) - centers(j,k))**2
                end do
                dist = sqrt(dist)
                
                if (dist < radius) then
                    ! Too close to an existing niche, update if better
                    if (fit(best_idx) > center_fits(j)) then
                        centers(j,:) = pop(best_idx,:)
                        center_fits(j) = fit(best_idx)
                    end if
                    is_new_niche = .false.
                    exit
                end if
            end do
            
            ! If it's far enough from all existing niches, add a new one
            if (is_new_niche) then
                do j = 1, num_nich
                    if (center_fits(j) < 0.0) then
                        centers(j,:) = pop(best_idx,:)
                        center_fits(j) = fit(best_idx)
                        exit
                    end if
                end do
            end if
        end if
    end subroutine update_niches

    ! Update niche centers based on current population
    subroutine update_niches_mahalanobis(pop, fit, centers, center_fits, num_nich, size, dims, radius, cov_valid, inv_cov, chi_threshold)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(size), intent(in) :: fit
        real, dimension(num_nich, dims), intent(inout) :: centers
        real, dimension(num_nich), intent(inout) :: center_fits
        real, dimension(num_nich, dims, dims), intent(in) :: inv_cov
        logical, dimension(num_nich), intent(in) :: cov_valid
        integer, intent(in) :: num_nich, size, dims
        real, intent(in) :: radius
        real, intent(in) :: chi_threshold
        integer :: i, j, k, best_idx
        real :: best_fit, dist
        logical :: is_new_niche


        ! First, find the best individual
        best_fit = -1.0
        best_idx = 1
        do i = 1, size
          if (fit(i) > best_fit) then
            best_fit = fit(i)
            best_idx = i
          end if
        end do
        
        ! Check if we need to add a new niche center
        if (center_fits(1) < 0.0) then
          ! First niche - just add it
          centers(1,:) = pop(best_idx,:)
          center_fits(1) = fit(best_idx)
        else
          ! Check if the best individual is close to any existing niche
          is_new_niche = .true.
          do j = 1, num_nich
            if (center_fits(j) < 0.0) exit  ! No more valid niches
            ! Use Mahalanobis distance preferentially if possible
            if (cov_valid(j)) then
              dist = mahalanobis_distance(pop(best_idx,:), centers(j,:),inv_cov(j,:,:),dims)
              if (dist**2 <= chi_threshold) then
                if (fit(best_idx) > center_fits(j)) then
                    centers(j,:) = pop(best_idx,:)
                    center_fits(j) = fit(best_idx)
                end if
                is_new_niche = .false.
                exit
              endif
            else 
              ! Or calculate Euclidean distance to this niche center
              dist = 0.0
              do k = 1, dims
                  dist = dist + (pop(best_idx,k) - centers(j,k))**2
              end do
              dist = sqrt(dist)
              if (dist < radius) then
                  ! Too close to an existing niche, update if better
                  if (fit(best_idx) > center_fits(j)) then
                      centers(j,:) = pop(best_idx,:)
                      center_fits(j) = fit(best_idx)
                  end if
                  is_new_niche = .false.
                  exit
              end if
            endif
          end do
          
          ! If it's far enough from all existing niches, add a new one
          if (is_new_niche) then
              do j = 1, num_nich
                  if (center_fits(j) < 0.0) then
                      centers(j,:) = pop(best_idx,:)
                      center_fits(j) = fit(best_idx)
                      exit
                  end if
              end do
          end if
        end if
    end subroutine update_niches_mahalanobis
    subroutine update_niches_improved(pop, fit, centers, center_fits, num_nich, size, dims, radius)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(size), intent(in) :: fit
        real, dimension(num_nich, dims), intent(inout) :: centers
        real, dimension(num_nich), intent(inout) :: center_fits
        integer, intent(in) :: num_nich, size, dims
        real, intent(in) :: radius
        
        ! Local variables
        integer :: i, j, k, idx
        integer, dimension(size) :: sorted_indices
        real :: dist
        logical :: is_new_niche, niche_updated
        integer :: candidates_to_consider, valid_niches
        
        ! Sort individuals by fitness (descending order)
        call sort_by_fitness(fit, sorted_indices, size)
        
        ! Consider top candidates (but not more than available niche slots)
        ! Process at least the top 10% of population or top 5 individuals, whichever is larger
        candidates_to_consider = max(5, size / 10)
        candidates_to_consider = min(candidates_to_consider, size)
        
        ! Count current valid niches
        valid_niches = 0
        do i = 1, num_nich
            if (center_fits(i) > -999.0) valid_niches = valid_niches + 1
        end do
        
        ! Process top candidates in order of fitness
        do i = 1, candidates_to_consider
            idx = sorted_indices(i)
            is_new_niche = .true.
            niche_updated = .false.
            
            ! Check against all existing niches
            do j = 1, num_nich
                if (center_fits(j) < -998.0) cycle  ! Skip invalid niches
                
                ! Calculate Euclidean distance to this niche center
                dist = 0.0
                do k = 1, dims
                    dist = dist + (pop(idx,k) - centers(j,k))**2
                end do
                dist = sqrt(dist)
                
                ! If close to existing niche
                if (dist < radius) then
                    is_new_niche = .false.
                    
                    ! Update niche if this individual is better
                    if (fit(idx) > center_fits(j)) then
                        centers(j,:) = pop(idx,:)
                        center_fits(j) = fit(idx)
                        niche_updated = .true.
                    end if
                    exit  ! Don't assign to multiple niches
                end if
            end do
            
            ! If it's a new niche and we have space, add it
            if (is_new_niche .and. valid_niches < num_nich) then
                ! Find first empty slot
                do j = 1, num_nich
                    if (center_fits(j) < -998.0) then
                        centers(j,:) = pop(idx,:)
                        center_fits(j) = fit(idx)
                        valid_niches = valid_niches + 1
                        exit
                    end if
                end do
            end if
            
            ! If we've filled all niche slots, we can stop early
            if (valid_niches >= num_nich) exit
        end do
    
    end subroutine update_niches_improved

! Helper subroutine to sort individuals by fitness (descending)
    subroutine sort_by_fitness(fitness, indices, size)
        integer, intent(in) :: size
        real, dimension(size), intent(in) :: fitness
        integer, dimension(size), intent(out) :: indices
        
        integer :: i, j, temp_idx
        real :: temp_fit
        real, dimension(size) :: fitness_copy
        
        ! Initialize indices
        do i = 1, size
            indices(i) = i
        end do
        
        ! Copy fitness array for sorting
        fitness_copy = fitness
        
        ! Simple bubble sort (could be replaced with quicksort for larger populations)
        do i = 1, size - 1
            do j = i + 1, size
                if (fitness_copy(i) < fitness_copy(j)) then
                    ! Swap fitness values
                    temp_fit = fitness_copy(i)
                    fitness_copy(i) = fitness_copy(j)
                    fitness_copy(j) = temp_fit
                    
                    ! Swap corresponding indices
                    temp_idx = indices(i)
                    indices(i) = indices(j)
                    indices(j) = temp_idx
                end if
            end do
        end do
    end subroutine sort_by_fitness
    
    ! Assign each individual to a niche using Mahalanobis distance where available
    ! Mahalonobis distance: doi.org/10.1007/s13171-019-00164-5
    ! Outlier analysis: 10.1007/978-3-319-47578-3
    ! Perfect explanation: https://www.cfholbert.com/blog/outlier_mahalanobis_distance/
    subroutine assign_niches_mahalanobis(pop, centers, inv_cov, cov_valid, membership, &
                                        size, num_nich, dims, radius, chi_threshold)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(num_nich, dims), intent(in) :: centers
        real, dimension(num_nich, dims, dims), intent(in) :: inv_cov
        logical, dimension(num_nich), intent(in) :: cov_valid
        integer, dimension(size), intent(out) :: membership
        integer, intent(in) :: size, num_nich, dims
        real, intent(in) :: radius, chi_threshold
        integer :: i, j, k, closest_niche
        real :: dist, min_dist

        
        do i = 1, size
            membership(i) = 0  ! Default - no niche
            min_dist = huge(min_dist)
            closest_niche = 0
            
            do j = 1, num_nich
                if (centers(j,1) < -998.0) cycle  ! Skip invalid niches
                
                if (cov_valid(j)) then
                    ! Use Mahalanobis distance
                    dist = mahalanobis_distance(pop(i,:), centers(j,:), inv_cov(j,:,:), dims)
                    ! Rescale it so the accepted values run between 0,2
                    ! Not ideal, but... at least a bit comparable?
                    dist = sqrt(2.0) * dist / sqrt(chi_threshold)
                else
                    ! Fall back to Euclidean distance
                    dist = 0.0
                    do k = 1, dims
                        dist = dist + (pop(i,k) - centers(j,k))**2
                    end do
                    dist = sqrt(dist)
                end if
          
                if (dist < min_dist) then
                    min_dist = dist
                    closest_niche = j
                end if
            end do
        ! Assign to closest niche with unified threshold
            if (closest_niche > 0) then
                if (cov_valid(closest_niche)) then
                    ! For Mahalanobis, use chi-square threshold
                    if (min_dist**2 <= 2.0) then
                        membership(i) = closest_niche
                    end if
                else
                    ! For Euclidean, use radius threshold
                    if (min_dist <= 2.0) then  ! Since we rescaled by radius
                        membership(i) = closest_niche
                    end if
                end if
            end if
        end do
    end subroutine assign_niches_mahalanobis
    
    ! Legacy assignment function (kept for compatibility, testing)
    subroutine assign_niches(pop, centers, membership, size, num_nich, dims, radius)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(num_nich, dims), intent(in) :: centers
        integer, dimension(size), intent(out) :: membership
        integer, intent(in) :: size, num_nich, dims
        real, intent(in) :: radius
        integer :: i, j, k, closest_niche
        real :: dist, min_dist
        
        do i = 1, size
            membership(i) = 0  ! Default - no niche
            min_dist = huge(min_dist)
            closest_niche = 0
            
            do j = 1, num_nich
                if (centers(j,1) < -998.0) cycle  ! Skip invalid niches
                
                dist = 0.0
                do k = 1, dims
                    dist = dist + (pop(i,k) - centers(j,k))**2
                end do
                dist = sqrt(dist)
                
                if (dist < min_dist) then
                    min_dist = dist
                    closest_niche = j
                end if
            end do
            
            ! Assign to closest niche if within radius
            if (closest_niche > 0 .and. min_dist <= 2.0 * radius) then
                membership(i) = closest_niche
            end if
        end do
    end subroutine assign_niches
    
    ! Create a new generation using selection, crossover, and mutation
    subroutine create_new_generation(pop, fit, new_pop, membership, size, dims, &
                                    num_nich, cross_prob, mut_prob, lower, upper)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(size), intent(in) :: fit
        real, dimension(size, dims), intent(out) :: new_pop
        integer, dimension(size), intent(in) :: membership
        integer, intent(in) :: size, dims, num_nich
        real, intent(in) :: cross_prob, mut_prob, lower, upper
        integer :: i, parent1, parent2
        real :: r
        real, dimension(size) :: shared_fitness
        
        ! Apply fitness sharing within niches
        call apply_fitness_sharing(pop, fit, membership, shared_fitness, size, dims, num_nich)
        
        ! Create new population
        do i = 1, size, 2
            ! Select two parents using tournament selection with shared fitness
            call tournament_select(shared_fitness, size, parent1)
            call tournament_select(shared_fitness, size, parent2)
            
            ! Crossover
            call random_number(r)
            if (r < cross_prob) then
                call crossover(pop(parent1,:), pop(parent2,:), new_pop(i,:), &
                              new_pop(min(i+1,size),:), dims, lower, upper)
            else
                new_pop(i,:) = pop(parent1,:)
                if (i+1 <= size) new_pop(i+1,:) = pop(parent2,:)
            end if
            
            ! Mutation
            call mutate(new_pop(i,:), dims, mut_prob, lower, upper)
            if (i+1 <= size) call mutate(new_pop(i+1,:), dims, mut_prob, lower, upper)
        end do
    end subroutine create_new_generation
    
    ! Apply fitness sharing to promote diversity within niches
    subroutine apply_fitness_sharing(pop, fit, membership, shared_fit, size, dims, num_nich)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(size), intent(in) :: fit
        integer, dimension(size), intent(in) :: membership
        real, dimension(size), intent(out) :: shared_fit
        integer, intent(in) :: size, dims, num_nich
        integer :: i, j, k, l, niche_count
        real :: dist, share_sum
        
        ! Initialize shared fitness to original fitness
        shared_fit = fit
        
        ! Apply sharing within each niche
        do i = 1, num_nich
            ! Count individuals in this niche
            niche_count = 0
            do j = 1, size
                if (membership(j) == i) niche_count = niche_count + 1
            end do
            
            if (niche_count > 1) then
                ! Apply sharing only if more than one individual in niche
                do j = 1, size
                    if (membership(j) == i) then
                      share_sum = 1.0  ! Start with self
                      
                      ! Add sharing component from other individuals in same niche
                      do k = 1, size
                          if (k /= j .and. membership(k) == i) then
                              dist = 0.0
                              do l = 1, dims
                                  dist = dist + (pop(j,l) - pop(k,l))**2
                              end do
                              dist = sqrt(dist)
                              
                              ! Triangular sharing function
                              if (dist < 0.1) then
                                  share_sum = share_sum + (1.0 - dist/0.1)
                              end if
                          end if
                      end do
                      
                      ! Apply sharing factor
                      shared_fit(j) = fit(j) / share_sum
                    end if
                end do
            end if
        end do
    end subroutine apply_fitness_sharing
    
    ! Tournament selection
    subroutine tournament_select(fit, size, selected)
        real, dimension(size), intent(in) :: fit
        integer, intent(in) :: size
        integer, intent(out) :: selected
        integer :: cand1, cand2
        real :: r
        
        call random_number(r)
        cand1 = int(r * size) + 1
        
        call random_number(r)
        cand2 = int(r * size) + 1
        
        if (fit(cand1) > fit(cand2)) then
            selected = cand1
        else
            selected = cand2
        end if
    end subroutine tournament_select
    
    ! Crossover operator (SBX - Simulated Binary Crossover)
    ! Agrawal, Ram & Deb, Kalyanmoy & Agrawal, Ram. (2000). Simulated Binary Crossover for Continuous Search Space. Complex Systems. 9. 
    ! Similar implementation to that of pymoo from MSU Coinlab: https://github.com/msu-coinlab/pymoo/pymoo/operators/crosssover/sbx
    ! DOI: 10.1109/ACCESS.2020.2990567 !
    subroutine crossover(parent1, parent2, child1, child2, dims, lower, upper)
        real, dimension(dims), intent(in) :: parent1, parent2
        real, intent(in) :: lower, upper
        real, dimension(dims), intent(out) :: child1, child2
        integer, intent(in) :: dims
        real :: beta, u
        integer :: i
        real, parameter :: eta_c = 20.0  ! Distribution index
        
        do i = 1, dims
            call random_number(u)
            if (u <= 0.5) then
                beta = (2.0 * u)**(1.0/(eta_c+1.0))
            else
                beta = (1.0/(2.0*(1.0-u)))**(1.0/(eta_c+1.0))
            end if
            
            child1(i) = 0.5 * ((1.0+beta)*parent1(i) + (1.0-beta)*parent2(i))
            child2(i) = 0.5 * ((1.0-beta)*parent1(i) + (1.0+beta)*parent2(i))
            ! Bound checking
            child1(i) = max(lower, min(upper, child1(i)))
            child2(i) = max(lower, min(upper, child2(i)))
        end do
    end subroutine crossover
    
    ! Mutation operator
    ! Highly-disruptive polynomial – 10.1016/j.ejor.2006.06.042, 10.1515/jisys-2018-0331
    subroutine mutate(indiv, dims, prob, lower, upper)
        integer, intent(in) :: dims
        real, dimension(dims), intent(inout) :: indiv
        real, intent(in) :: prob, lower, upper
        integer :: i
        real :: r, delta, delta_1, delta_2, delta_max
        real, parameter :: eta_m = 20.0  ! Distribution index
        
        do i = 1, dims
            call random_number(r)
            if (r < prob) then

                delta_max = (upper - lower)
                delta_1 = (indiv(i) - lower) / delta_max
                delta_2 = (upper - indiv(i)) / delta_max
                call random_number(r) 
                
                if (r < 0.5) then
                    delta = (2.0*r + (1.0 - 2.0*r)*&
                           &(1.0-delta_1)**(eta_m+1.0))**((1.0/(eta_m+1.0)))-1
                else
                    delta = 1.0 - (2.0*(1.0-r)+2.0*(r-0.5)*&
                           &(1.0-delta_2)**(eta_m+1.0))**(1.0/(eta_m+1.0))
                end if  
               
                ! Apply mutation
                indiv(i) = indiv(i) + delta * delta_max
                
                ! Bound checking
                indiv(i) = max(lower, min(upper, indiv(i)))
            end if
        end do
    end subroutine mutate
    
end program CumulativeMultiNichingGA
