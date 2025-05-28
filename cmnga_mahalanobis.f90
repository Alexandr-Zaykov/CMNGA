program CumulativeMultiNichingGA
    implicit none
    
    ! External BLAS/LAPACK functions
    external :: sgemv, sdot, sger, sgetrf, sgetri, sgecon, spotrf, spotri
    real :: sdot
    
    ! Parameters
    ! There are further distribution parameters for crossover and mutation "eta_c" – we use 20 for both for now.
    integer, parameter :: pop_size = 40       ! Population size
    integer, parameter :: max_gen = 200          ! Maximum generations
    integer, parameter :: num_vars = 2           ! Number of variables in optimization (up to 12D)
    integer, parameter :: num_niches = 15        ! Number of niches to maintain
    real, parameter :: mutate_rate = 0.30         ! Mutation rate
    real, parameter :: cross_rate = 0.05           ! Crossover rate
    real, parameter :: radius_factor = 0.05       ! Niche radius factor
    real, parameter :: min_bound = -4.0           ! Lower bound for variables
    real, parameter :: max_bound = 4.0            ! Upper bound for variables
    real, parameter :: min_niche_size = 3         ! Minimum individuals per niche for covariance
    real, parameter :: reg_param = 1e-6           ! Regularization parameter for covariance
    
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
        call update_niches(population, fitness, niche_centers, niche_fitness, &
                          num_niches, pop_size, num_vars, niche_radius)
        
        ! Assign individuals to niches
        call assign_niches(population, niche_centers, niche_membership, &
                          pop_size, num_niches, num_vars, niche_radius)
        
        ! Update covariance matrices for each niche
        call update_niche_covariances(population, niche_membership, niche_centers, &
                                     niche_covariance, niche_inv_cov, cov_valid, &
                                     pop_size, num_niches, num_vars)
        
        ! Re-assign individuals using Mahalanobis distance where possible
        call assign_niches_mahalanobis(population, niche_centers, niche_inv_cov, &
                                      cov_valid, niche_membership, pop_size, &
                                      num_niches, num_vars, niche_radius)
        
        ! Create new population
        call create_new_generation(population, fitness, new_population, niche_membership, &
                                  pop_size, num_vars, num_niches, cross_rate, mutate_rate, &
                                  min_bound, max_bound)
        
        ! Replace old population with new one
        population = new_population
        
        ! Print progress
        if (mod(gen, 10) == 0) then
            print *, "Generation ", gen
            print *, "Best niches found:"
            do i = 1, num_niches
                if (niche_fitness(i) > -999.0) then
                    print *, "Niche ", i, ": Position = ", niche_centers(i,:), &
                           " Fitness = ", niche_fitness(i), " Cov_valid = ", cov_valid(i)
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
            print *, "Niche ", i, ": Position = ", niche_centers(i,:), &
                   " Value = ", niche_fitness(i) - 1.0, " Cov_valid = ", cov_valid(i)
        end if
    end do
    call CPU_TIME(time_end)
    print *, "Time taken: ", time_end - time_start, " s"
    print *, "Anticipated Results: +-0.31,0.93,1.55,2.17,2.79,3.41"
    print *, "Number of evaluations: ", num_eval
    
contains

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
       
        ! ADD THE SEARCH FUNCTION HERE
        do i = 1, size
            evalfs = evalfs + 1
            ! Example multi-modal function: Sum of sines
            if (dims == 1) then
                fit(i) = sin(5.0 * pop(i,1))**6 * exp(-pop(i,1)**2)
                ! Results: +-0.31,0.93,1.55,2.17,2.79,3.41
            else
                ! 2D or higher: Combination of sines (only use first 2 dimensions for test)
                fit(i) = sin(5.0 * pop(i,1))**2 * sin(5.0 * pop(i,2))**2 * &
                         exp(-(pop(i,1)**2 + pop(i,2)**2)/2)
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
        
        ! dist = sqrt(diff^T * temp) (using BLAS SDOT)
        dist = sqrt(max(0.0, sdot(dims, diff, 1, temp, 1)))
    end function mahalanobis_distance
    
    ! Update covariance matrices for each niche using LAPACK
    subroutine update_niche_covariances(pop, membership, centers, cov_matrices, &
                                       inv_cov_matrices, cov_valid, size, num_nich, dims)
        real, dimension(size, dims), intent(in) :: pop
        integer, dimension(size), intent(in) :: membership
        real, dimension(num_nich, dims), intent(in) :: centers
        real, dimension(num_nich, dims, dims), intent(out) :: cov_matrices
        real, dimension(num_nich, dims, dims), intent(out) :: inv_cov_matrices
        logical, dimension(num_nich), intent(out) :: cov_valid
        integer, intent(in) :: size, num_nich, dims
        
        integer :: niche_id, i, j, k, count, info
        real, dimension(dims) :: mean_vec, diff
        real, dimension(dims, dims) :: cov_mat, work_mat
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
            if (count >= int(min_niche_size)) then
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
                
                ! Store covariance matrix
                cov_matrices(niche_id, :, :) = cov_mat
                
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
    
    ! Assign each individual to a niche using Mahalanobis distance where available
    subroutine assign_niches_mahalanobis(pop, centers, inv_cov, cov_valid, membership, &
                                        size, num_nich, dims, radius)
        real, dimension(size, dims), intent(in) :: pop
        real, dimension(num_nich, dims), intent(in) :: centers
        real, dimension(num_nich, dims, dims), intent(in) :: inv_cov
        logical, dimension(num_nich), intent(in) :: cov_valid
        integer, dimension(size), intent(out) :: membership
        integer, intent(in) :: size, num_nich, dims
        real, intent(in) :: radius
        integer :: i, j, k, closest_niche
        real :: dist, min_dist
        real :: chi_threshold
        
        do i = 1, size
            membership(i) = 0  ! Default - no niche
            min_dist = huge(min_dist)
            closest_niche = 0
            
            do j = 1, num_nich
                if (centers(j,1) < -998.0) cycle  ! Skip invalid niches
                
                if (cov_valid(j)) then
                    ! Use Mahalanobis distance
                    dist = mahalanobis_distance(pop(i,:), centers(j,:), inv_cov(j,:,:), dims)
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
            
            ! Assign to closest niche if within radius (use adaptive radius for Mahalanobis)
            if (closest_niche > 0) then
                if (cov_valid(closest_niche)) then
                    ! For Mahalanobis distance, use chi-square threshold
                    ! For dims degrees of freedom, 95% confidence interval
                    ! This needs checking
                    if (dims <= 1) then
                        chi_threshold = 3.84  ! chi2(0.95, 1)
                    else if (dims <= 2) then
                        chi_threshold = 5.99  ! chi2(0.95, 2)
                    else if (dims <= 3) then
                        chi_threshold = 7.81  ! chi2(0.95, 3)
                    else if (dims <= 4) then
                        chi_threshold = 9.49  ! chi2(0.95, 4)
                    else if (dims <= 5) then
                        chi_threshold = 11.07 ! chi2(0.95, 5)
                    else if (dims <= 6) then
                        chi_threshold = 12.59 ! chi2(0.95, 6)
                    else if (dims <= 8) then
                        chi_threshold = 15.51 ! chi2(0.95, 8)
                    else if (dims <= 10) then
                        chi_threshold = 18.31 ! chi2(0.95, 10)
                    else
                        chi_threshold = 21.03 ! chi2(0.95, 12)
                    end if
                    
                    if (min_dist**2 <= chi_threshold) then  ! Square for comparison
                        membership(i) = closest_niche
                    end if
                else
                    ! Use original Euclidean threshold
                    if (min_dist <= 2.0 * radius) then
                        membership(i) = closest_niche
                    end if
                end if
            end if
        end do
    end subroutine assign_niches_mahalanobis
    
    ! Legacy assignment function (kept for compatibility)
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
