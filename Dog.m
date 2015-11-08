classdef Dog < handle
    %   A boid, hanlde subclass means reference type
    
    properties
        position = zeros(1,2)
        velocity = zeros(1,2)
        sheepMass = zeros(1,4)% number of sheep in each quadrant
        deltaVelocity = zeros(1,2)
        sheepMassHistory = [] %history of sheepMass over a sim
        deltaVelocityHistory = [] %history of delta velocity over a sim
        positionHistory = []
        sheepPosHistory = []
        velocityHistory = []
    end
    
    methods
        % Automatically store history of variables.
        function obj = set.sheepMass(obj, newValue)
            obj.sheepMassHistory = [obj.sheepMassHistory; 
                                                newValue];
        end
        function obj = set.deltaVelocity(obj, newValue)
            obj.deltaVelocityHistory = [obj.deltaVelocityHistory; 
                                                        newValue];
        end
        function obj = set.position(obj, newValue)
            obj.positionHistory = [obj.positionHistory; 
                                              newValue];
        end
        function obj = set.velocity(obj, newValue)
            obj.velocityHistory = [obj.velocityHistory; 
                                             newValue];
        end
    end
    
end
