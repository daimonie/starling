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
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy, speedsquared) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i < j) then
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distancesquared < sensitivities(i)**2 ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        newvelocities(j,:) = newvelocities(j,:) + velocities(i,:)
                    end if
                end if
            end do
        end do
        !$omp end parallel do
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i, speedsquared, alpha, beta, gamma, rotationz, rotationx, &
        !$omp rotationy, rotation) &
        !$omp firstprivate(eta,pi, velocities) &
        !$omp shared(newvelocities)
        do i = 1, number 
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
    
    subroutine bowlhabitat(number, positions, velocities, sensitivities, tau, eta, newpositions, &
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
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i < j) then
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distancesquared < sensitivities(i)**2 ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        newvelocities(j,:) = newvelocities(j,:) + velocities(i,:)
                    end if
                end if
            end do
            !Now, distance from centre
            distancesquared = dot_product(positions(i,:), positions(i,:))
            if (distancesquared**0.5 + sensitivities(i) > habitatsize) then
                newvelocities(i,:) = newvelocities(i,:) - habitatstrength * positions(i,:)/distancesquared**0.5
            end if  
        end do
        !$omp end parallel do
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i, speedsquared, alpha, beta, gamma, rotationz, rotationx, &
        !$omp rotationy, rotation) &
        !$omp firstprivate(eta,pi, velocities) &
        !$omp shared(newvelocities)
        do i = 1, number 
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
    end subroutine bowlhabitat
    subroutine interactionbowlhabitat(number, positions, velocities, sensitivities, tau, eta, i0, i1, &
        i2, i3, i4, i5, newpositions, newvelocities, habitatsize, habitatstrength)
        implicit none 
        !Pretty much the same function but we simulate a spherical habitat by using virtual fish
        integer, intent(in) :: number
        double precision, intent(in):: tau, eta, habitatSize, habitatstrength, &
            i0, i1, i2, i3, i4, i5
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(number) :: sensitivities

        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared, force, x
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared, force, x) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength, &
        !$omp i0, i1, i2, i3, i4, i5) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i < j) then
                    force = 0.00_8
                        
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    x = distancesquared**0.5
                    if (x < sensitivities(i) ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        newvelocities(j,:) = newvelocities(j,:) + velocities(i,:)
                        
                        !Let us add an unnecessarily complex calculation for some sort of force.  
                        !The idea is that it repulses at close distance, attracts at distances near sensitivity_i 
                    else if( x < i3*3.0) then  
!                         force = i0 + i1 / (i2 + x) + i4 * (x - i3) / (1.0 + exp( (x-sensitivities(i))/(2.0*i5**2))) 
                    !     y[i] = i0 + i1/( 1.0 + np.exp((thisx-i1)/(2*i2**2))) 
                    !     if thisx > i3:
                    !         y[i] += i4*(thisx-i3)*1.0 / ( 1.0 + np.exp( (thisx-s)/(2*i5**2)))
                        force = i0 + i1 / (1.0 + exp( (x-i1)/(2*i2**2)))
                        if (x > i3) then
                            force = force + i4 * (x-i3) * 1.0 / (1.0 + exp( (x-sensitivities(i))/(2*i5*2)))
                        end if
                    end if
                    newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - positions(j,:))/ x * force * tau
                    newvelocities(j,:) = newvelocities(j,:) - (positions(i,:) - positions(j,:))/ x * force * tau
                end if
            end do
            !Now, distance from centre
            distancesquared = dot_product(positions(i,:), positions(i,:))
            if (distancesquared**0.5 + sensitivities(i) > habitatsize) then
                newvelocities(i,:) = newvelocities(i,:) - habitatstrength * positions(i,:)/distancesquared**0.5
            end if 
        end do
        !$omp end parallel do
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i, speedsquared, alpha, beta, gamma, rotationz, rotationx, &
        !$omp rotationy, rotation) &
        !$omp firstprivate(eta,pi, velocities) &
        !$omp shared(newvelocities)
        do i = 1, number 
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
    end subroutine interactionbowlhabitat
    subroutine interactionboxhabitat(number, positions, velocities, sensitivities, tau, eta, i0, i1, &
        i2, i3, i4, i5, newpositions, newvelocities, habitatsize, habitatstrength)
        implicit none 
        !Pretty much the same function but we simulate a spherical habitat by using virtual fish
        integer, intent(in) :: number
        double precision, intent(in):: tau, eta, habitatSize, habitatstrength, &
            i0, i1, i2, i3, i4, i5
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(number) :: sensitivities

        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared, force, x
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared, force, x) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength, &
        !$omp i0, i1, i2, i3, i4, i5) &
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i < j) then
                    force = 0.00_8
                        
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    x = distancesquared**0.5
                    if (x < sensitivities(i) ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        newvelocities(j,:) = newvelocities(j,:) + velocities(i,:)
                        
                        !Let us add an unnecessarily complex calculation for some sort of force.  
                        !The idea is that it repulses at close distance, attracts at distances near sensitivity_i 
                    else if( x < i3*3.0) then  
!                         force = i0 + i1 / (i2 + x) + i4 * (x - i3) / (1.0 + exp( (x-sensitivities(i))/(2.0*i5**2))) 
                    !     y[i] = i0 + i1/( 1.0 + np.exp((thisx-i1)/(2*i2**2))) 
                    !     if thisx > i3:
                    !         y[i] += i4*(thisx-i3)*1.0 / ( 1.0 + np.exp( (thisx-s)/(2*i5**2)))
                        force = i0 + i1 / (1.0 + exp( (x-i1)/(2*i2**2)))
                        if (x > i3) then
                            force = force + i4 * (x-i3) * 1.0 / (1.0 + exp( (x-sensitivities(i))/(2*i5*2)))
                        end if
                    end if
                    newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - positions(j,:))/ x * force * tau
                    newvelocities(j,:) = newvelocities(j,:) - (positions(i,:) - positions(j,:))/ x * force * tau
                end if
            end do
            !Now, distance from centre
            distancesquared = dot_product(positions(i,:), positions(i,:))
            if (distancesquared**0.5 + sensitivities(i) > habitatsize) then
                newvelocities(i,:) = newvelocities(i,:) - habitatstrength * positions(i,:)/distancesquared**0.5
            end if 
        end do
        !$omp end parallel do
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i, speedsquared, alpha, beta, gamma, rotationz, rotationx, &
        !$omp rotationy, rotation) &
        !$omp firstprivate(eta,pi, velocities) &
        !$omp shared(newvelocities)
        do i = 1, number 
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
    end subroutine interactionboxhabitat
    
    subroutine sharkbowl(number, positions, velocities, sensitivities, tau, eta, i0, i1, &
        i2, i3, i4, i5, newpositions, newvelocities, habitatsize, habitatstrength, &
        predatorsense, predatorstrength, predatorlocation, predatornumber)
        implicit none 
        !Pretty much the same function but we simulate a spherical habitat by using virtual fish
        integer, intent(in) :: number, predatornumber
        double precision, intent(in):: tau, eta, habitatSize, habitatstrength, &
            i0, i1, i2, i3, i4, i5
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(predatornumber, 3) :: predatorLocation
        double precision, intent(in), dimension(number) :: sensitivities
        double precision, intent(in), dimension(predatornumber) :: predatorSense, predatorStrength

        
        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j, p
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared, force, x
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared, force,x,p) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength, &
        !$omp i0, i1, i2, i3, i4, i5, &
        !$omp predatorsense, predatorstrength, predatorlocation, predatornumber)&
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i < j) then
                    force = 0.00_8
                        
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    x = distancesquared**0.5
                    if (x < sensitivities(i) ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        newvelocities(j,:) = newvelocities(j,:) + velocities(i,:)
                        
                        !Let us add an unnecessarily complex calculation for some sort of force.  
                        !The idea is that it repulses at close distance, attracts at distances near sensitivity_i 
                    else if( x < i3*3.0) then  
!                         force = i0 + i1 / (i2 + x) + i4 * (x - i3) / (1.0 + exp( (x-sensitivities(i))/(2.0*i5**2))) 
                    !     y[i] = i0 + i1/( 1.0 + np.exp((thisx-i1)/(2*i2**2))) 
                    !     if thisx > i3:
                    !         y[i] += i4*(thisx-i3)*1.0 / ( 1.0 + np.exp( (thisx-s)/(2*i5**2)))
                        force = i0 + i1 / (1.0 + exp( (x-i1)/(2*i2**2)))
                        if (x > i3) then
                            force = force + i4 * (x-i3) * 1.0 / (1.0 + exp( (x-sensitivities(i))/(2*i5*2)))
                        end if
                    end if
                    newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - positions(j,:))/ x * force * tau
                    newvelocities(j,:) = newvelocities(j,:) - (positions(i,:) - positions(j,:))/ x * force * tau
                end if
            end do
            do p = 1, predatornumber
                distancesquared = dot_product( positions(i,:) - predatorlocation(p,:), positions(i,:) - predatorlocation(p,:))
                x = distancesquared**0.5
                if (x < predatorsense(i)) then 
                    newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - predatorlocation(p,:))/ x * predatorstrength(p) 
                end if
            end do
            !Now, distance from centre
            distancesquared = dot_product(positions(i,:), positions(i,:))
            if (distancesquared**0.5 + sensitivities(i) > habitatsize) then
                newvelocities(i,:) = newvelocities(i,:) - habitatstrength * positions(i,:)/distancesquared**0.5
            end if 
        end do
        !$omp end parallel do
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i, speedsquared, alpha, beta, gamma, rotationz, rotationx, &
        !$omp rotationy, rotation) &
        !$omp firstprivate(eta,pi, velocities) &
        !$omp shared(newvelocities)
        do i = 1, number 
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
    end subroutine sharkbowl  
    
    subroutine sharkbox(number, positions, velocities, sensitivities, tau, eta, i0, i1, &
        i2, i3, i4, i5, newpositions, newvelocities, habitatsize, habitatstrength, &
        predatorsense, predatorstrength, predatorlocation, predatornumber)
        implicit none 
        !Pretty much the same function but we simulate a spherical habitat by using virtual fish
        integer, intent(in) :: number, predatornumber
        double precision, intent(in):: tau, eta, habitatSize, habitatstrength, &
            i0, i1, i2, i3, i4, i5
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(predatornumber, 3) :: predatorLocation
        double precision, intent(in), dimension(number) :: sensitivities
        double precision, intent(in), dimension(predatornumber) :: predatorSense, predatorStrength

        
        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities

        integer :: i, j, p
        double precision :: distancesquared, pi, alpha, beta, gamma, speedsquared, force, x
        
        double precision, dimension(3,3) :: rotation, rotationx, rotationy, rotationz
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads  
        
        pi = 4.0*atan(1.004)
        
        newpositions = positions
        newvelocities = velocities

        distancesquared = 0.00
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared,alpha,beta,gamma,rotation,rotationz,rotationx,rotationy,speedsquared, force,x,p) &
        !$omp firstprivate(number,positions,velocities,tau,eta,sensitivities,pi, habitatsize, habitatstrength, &
        !$omp i0, i1, i2, i3, i4, i5, &
        !$omp predatorsense, predatorstrength, predatorlocation, predatornumber)&
        !$omp shared(newpositions,newvelocities)
        do i = 1, number
            newpositions(i,:) = positions(i,:) + tau * velocities(i, :)
            newvelocities(i,:) = velocities(i,:)
            do j = 1, number
                if ( i < j) then 
                    force = 0.00_8
                        
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    x = distancesquared**0.5
                    if (x < sensitivities(i) ) then
                        newvelocities(i,:) = newvelocities(i,:) + velocities(j,:)
                        newvelocities(j,:) = newvelocities(j,:) + velocities(i,:)
                        
                        !Let us add an unnecessarily complex calculation for some sort of force.  
                        !The idea is that it repulses at close distance, attracts at distances near sensitivity_i 
                    else if( x < i3*3.0) then   
                        force = i0 + i1 / (1.0 + exp( (x-i1)/(2*i2**2)))
                        if (x > i3) then
                            force = force + i4 * (x-i3) * 1.0 / (1.0 + exp( (x-sensitivities(i))/(2*i5*2)))
                        end if
                    end if
                    newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - positions(j,:))/ x * force * tau
                    newvelocities(j,:) = newvelocities(j,:) - (positions(i,:) - positions(j,:))/ x * force * tau
                end if
                !Attraction to centre, like a spring
                x = ( dot_product( positions(i,:),positions(i,:)))**0.5
                force = -habitatstrength / habitatsize * x
                newvelocities(i,:) = newvelocities(i,:) + positions(i,:)/ x * force * tau
            end do
            do p = 1, predatornumber
                distancesquared = dot_product( positions(i,:) - predatorlocation(p,:), positions(i,:) - predatorlocation(p,:))
                x = distancesquared**0.5
                if (x < predatorsense(i)) then 
                    newvelocities(i,:) = newvelocities(i,:) + (positions(i,:) - predatorlocation(p,:))/ x * predatorstrength(p) 
                end if
            end do 
        end do
        !$omp end parallel do
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i, speedsquared, alpha, beta, gamma, rotationz, rotationx, &
        !$omp rotationy, rotation) &
        !$omp firstprivate(eta,pi, velocities) &
        !$omp shared(newvelocities)
        do i = 1, number 
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
    end subroutine sharkbox 
    
    subroutine rotatepoints(number, points, axis, theta, returnpoints)
        implicit none
        
        integer, intent(in) :: number
        
        double precision, intent(in), dimension(number, 3) :: points 
        double precision, intent(in), dimension(3) :: axis
        double precision, intent(in) :: theta
        
        double precision, intent(out), dimension(number, 3) :: returnpoints
      
        double precision, dimension(3) :: centre
        double precision, dimension(3,3) :: rotation 
        
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads, i
        
        centre = sum(points, dim=1) / number 
         
        rotation(1,1) = cos(theta) + axis(1)**2 * ( 1 - cos(theta))
        rotation(1,2) = axis(1) * axis(2) * (1 - cos(theta)) - axis(3) * sin(theta)
        rotation(1,3) = axis(1) * axis(3) * (1 - cos(theta)) + axis(2) * sin(theta)
        
        
        rotation(2,1) = axis(2) * axis(1) * (1 - cos(theta)) + axis(3) * sin(theta)
        rotation(2,2) = cos(theta) + axis(2)**2 * (1 - cos(theta))
        rotation(2,3) = axis(2) * axis(3) * (1 - cos(theta)) - axis(1) * sin(theta)
        
        rotation(3,1) = axis(3) * axis(1) * (1 - cos(theta)) - axis(2) * sin(theta)
        rotation(3,2) = axis(3) * axis(2) * (1 - cos(theta)) + axis(1) * sin(theta)
        rotation(3,3) = cos(theta) + axis(3)**2 * (1 - cos(theta))
        print *, "Rotation ends up at ", sum( matmul(rotation, transpose(rotation))) 
        
        
        call omp_set_num_threads(omp_get_max_threads());  
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i) &
        !$omp firstprivate(centre, rotation, points) &
        !$omp shared(returnpoints)
        do i = 1, number
            returnpoints(i,:) = centre + matmul( rotation, points(i,:) - centre)
        end do 
        !$omp end parallel do
    end subroutine rotatepoints

    subroutine simplerotation(number, positions, orientations, viscosity, cutOff, tau, Tboundary1, Tboundary2, Tboundary3, &
        Talign, neworientations)!,A)
        implicit none
        
        integer, intent(in) :: number
	double precision, intent(in) :: viscosity, cutOff, tau, Tboundary1, Tboundary2, Tboundary3, Talign	
        
        double precision, intent(in), dimension(number, 3) :: positions, orientations 

        double precision, intent(out), dimension(number, 3) :: neworientations
	!double precision, intent(out), dimension(number, 3) :: A
     
        integer :: omp_get_thread_num, id, threadNum, omp_get_max_threads, i, j
	double precision :: distancesquared, factor
	double precision, dimension(number) :: numNeighbors
	double precision, dimension(number) :: angleCM, localCMmagnitude, orientationsMagnitude
	double precision, dimension(number, 3) :: torque, torqueBoundary, localCM, normOrient, torqueAlign

	do i=1,number
		normOrient(i,:) = orientations(i,:)/sqrt(dot_product(orientations(i,:),orientations(i,:)))
	end do

	torqueAlign(:,:) = 0.00
	distancesquared = 0.00
	numNeighbors(:) = 0
	localCM(:,:) = 0.00
	call omp_set_num_threads(omp_get_max_threads());
        !$omp parallel do &
        !$omp default(none) & 
        !$omp private(i,j,distancesquared) &
        !$omp firstprivate(number,positions,cutOff) &
        !$omp shared(torqueAlign,localCM,numNeighbors,normOrient)
        do i = 1, number
            do j = 1, number
                if (i < j) then
                    distancesquared = dot_product( positions(i,:) - positions(j,:), positions(i,:) - positions(j,:))
                    if (distancesquared < cutOff**2 ) then
                        torqueAlign(i,:) = torqueAlign(i,:) + normOrient(i,:) - normOrient(j,:)
                        torqueAlign(j,:) = torqueAlign(j,:) + normOrient(j,:) - normOrient(i,:)
			localCM(i,:) = localCM(i,:) + positions(j,:)
			localCM(j,:) = localCM(j,:) + positions(i,:)
			numNeighbors(i) = numNeighbors(i) + 1
			numNeighbors(j) = numNeighbors(j) + 1
                    end if
                end if
            end do
        end do	
        !$omp end parallel do

	factor = 0.00
        torqueAlign = torqueAlign*(-1.0)*Talign
	torqueBoundary(:,:) = 0.00
	!For each paricle, calculate the angle between the orientation of the particle and the vector pointing to the center of mass of its domain of vision
	!Also calculate torqueBoundary, the torque that tries to align the orientation of each particle in the direction of the flock bulk. Three parameters can be varied
	angleCM(:) = 0.00
	do i = 1, number
		if (numNeighbors(i) > 0) then
			localCM(i,:) = localCM(i,:) / numNeighbors(i)
			localCMmagnitude(i) = sqrt(dot_product(localCM(i,:), localCM(i,:)))
			orientationsMagnitude(i) = sqrt(dot_product(normOrient(i,:), normOrient(i,:)))
			angleCM(i) = acos( dot_product(normOrient(i,:), localCM(i,:)) / (localCMmagnitude(i)*orientationsMagnitude(i)) )
			factor = Tboundary1 * angleCM(i) + Tboundary2 * localCMmagnitude(i) + Tboundary3 * angleCM(i) * localCMmagnitude(i)
			torqueBoundary(i, :) = (localCM(i,:) / localCMmagnitude(i)) * factor
		else if (numNeighbors(i) == 0) then
			torqueBoundary(i, :) = [0.0,0.0,0.0]
		end if
	end do
	
	torque = torqueAlign + torqueBoundary

	neworientations = normOrient + tau * torque / viscosity
	do i = 1, number
		neworientations(i,:) = neworientations(i,:) / (dot_product(neworientations(i,:),neworientations(i,:)))**0.5
	end do
	!A(:,:) =0.0
    end subroutine simplerotation
end module ethology
