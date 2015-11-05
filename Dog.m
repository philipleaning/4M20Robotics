classdef Dog < handle
    %   A boid, hanlde subclass means reference type
    
    properties
        position = zeros(1,2)
        velocity = zeros(1,2)
        sheepMass = zeros(1,4)% number of sheep in each quadrant
        deltaVelocity = zeros(1,2)
        sheepMassHistory = zeros(20000,4) %history of sheepMass over a sim
        deltaVelocityHistory = zeros(20000,2) %history of delta velocity over a sim
        positionHistory = zeros(20000,2)
        sheepPosHistory = zeros(20000,20)
        velocityHistory = zeros(20000,2)
    end
    
    methods
        
    end
    
end
