!f90
!f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -m -c fortranMolDy mol-dy.f90; 
!@Anton: ethology is the science of animal behaviour. It's a sensible name.
module ethology
    contains
    subroutine simplebehaviour(number, positions, velocities, sensitivities, tau, eta, newpositions, newvelocities)
        implicit none

        integer, intent(in) :: number
        double precision, intent(in):: tau, eta
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(number) :: sensitivities

        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j
        double precision :: distanceSquared, pi, alpha, beta, gamma
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads;  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distanceSquared = 0.00
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distanceSquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i /= j) then
                    distanceSquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distanceSquared < sensitivities(i)**2 ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                    end if
                end if
            end do
            
            newvelocities(i,:) = newvelocities(i,:) / (dot_product(newvelocities(i,:), newvelocities(i,:)) )**0.5
             
            alpha = eta * pi - 2 * eta * pi * rand() 
            beta  = eta * pi - 2 * eta * pi * rand() 
            gamma = eta * pi - 2 * eta * pi * rand() 
            
            rotationz = reshape([1.0_8,0.0_8,0.0_8,0.0_8,cos(alpha),sin(alpha),0.0_8,-sin(alpha),cos(alpha)],shape(rotationz))
            rotationy = reshape([cos(beta),0.0_8,-sin(beta),0.0_8,1.0_8,0.0_8,sin(beta),0.0_8,cos(beta)], shape(rotationy))
            rotationx = reshape([cos(gamma),sin(gamma),0.0_8,-sin(gamma),cos(gamma),0.0_8,0.0_8,0.0_8,1.0_8],shape(rotationx))
            
            rotation = rotationz * rotationy * rotationx 
            newvelocities(i,:) =  matmul(rotation, newvelocities(i,:))
        end do
        !$omp end parallel do
    end subroutine simplebehaviour
end module ethology