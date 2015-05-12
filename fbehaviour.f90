!f90
!f2py --fcompiler=gfortran --f90flags="-fopenmp" -lgomp -m -c fortranMolDy mol-dy.f90; 
!@Anton: ethology is the science of animal behaviour. It's a sensible name.
module ethology
    contains
    subroutine simplebehaviour(number, positions, velocities, sensitivities, newpositions, newvelocities)
        implicit none
        integer, intent(in) :: number
        double precision, intent(in), dimension(number, 3) :: positions, velocities
        double precision, intent(in), dimension(number, 1) :: sensitivities
        double precision, intent(out), dimension(number, 3) :: newpositions, newvelocities
        
        newpositions = positions
        newvelocities = velocities
        
    end subroutine simplebehaviour
end module ethology