classdef Dog < handle
    %   A boid, hanlde subclass means reference type
    
    properties
        position = zeros(1,2)
        velocity = zeros(1,2)
        sheepMass = zeros(1,4)% number of sheep in each quadrant
        deltaVelocity = zeros(1,2)
        sheepMassHistory = [] %history of sheepMass over a sim
        deltaVelocityHistory = [] %history of delta velocity over a sim
    end
    
    methods
        
    end
    
end
