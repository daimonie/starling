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
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads;  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy, speedsquared) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i /= j) then
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distancesquared < sensitivities(i)**2 ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                    end if
                end if
            end do
             
            speedsquared = dot_product(newvelocities(i,:), newvelocities(i,:))
            if (speedsquared > 0) then
                newvelocities(i,:) = newvelocities(i,:) / (speedsquared)**0.5
            else
                newvelocities(i,:) = velocities(i,:)
            end if 
            
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
    
    subroutine simplehabitat(number, positions, velocities, sensitivities, tau, eta, newpositions, &
        newvelocities, habitatsize, habitatstrength)
        implicit none
        !Pretty much the same function but we simulate a spherical habitat by using virtual fish
        integer, intent(in) :: number
        double precision, intent(in):: tau, eta, habitatSize, habitatstrength
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(number) :: sensitivities

        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads;  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i /= j) then
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distancesquared < sensitivities(i)**2 ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                    end if
                end if
            end do
            !Now, distance from centre
            distancesquared = dot_product(positions(i,:), positions(i,:))
            if (distancesquared**0.5 + sensitivities(i) > habitatsize) then
                newvelocities(i,:) = newvelocities(i,:) - habitatstrength * positions(i,:)/distancesquared**0.5
            end if 
            speedsquared = dot_product(newvelocities(i,:), newvelocities(i,:))
            if (speedsquared > 0) then
                newvelocities(i,:) = newvelocities(i,:) / (speedsquared)**0.5
            else
                newvelocities(i,:) = velocities(i,:)
            end if
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
    end subroutine simplehabitat
    subroutine interactionhabitat(number, positions, velocities, sensitivities, tau, eta, i0, i1, &
        i2, i3, newpositions, newvelocities, habitatsize, habitatstrength)
        implicit none
        !Pretty much the same function but we simulate a spherical habitat by using virtual fish
        integer, intent(in) :: number
        double precision, intent(in):: tau, eta, habitatSize, habitatstrength, &
            i0, i1, i2, i3
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(number) :: sensitivities

        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared, force, x
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads;  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared, force, x) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength, &
        !$omp i0, i1, i2, i3) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i /= j) then
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distancesquared < sensitivities(i)**2 ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        
                        !Let us add an unnecessarily complex calculation for some sort of force.  
                        !The idea is that it repulses at close distance, attracts at distances near sensitivity_i
                        x = distanceSquared**0.5
                        force = 0.00_8
                        if (x < i0) then
                            force = i1
                        else if (x > i2) then
                            force = i3
                        end if 
                        
                        newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - positions(j,:))/ x * force * tau
                        
                    end if
                end if
            end do
            !Now, distance from centre
            distancesquared = dot_product(positions(i,:), positions(i,:))
            if (distancesquared**0.5 + sensitivities(i) > habitatsize) then
                newvelocities(i,:) = newvelocities(i,:) - habitatstrength * positions(i,:)/distancesquared**0.5
            end if
            speedsquared = dot_product(newvelocities(i,:), newvelocities(i,:))
            if (speedsquared > 0) then
                newvelocities(i,:) = newvelocities(i,:) / (speedsquared)**0.5
            else
                newvelocities(i,:) = velocities(i,:)
            end if
            
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
    end subroutine interactionhabitat
end module ethology