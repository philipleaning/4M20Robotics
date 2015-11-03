classdef Boid < handle
    %   A boid, hanlde subclass means reference type
    
    properties
        position = zeros(1,2)
        velocity = zeros(1,2)
        % If afraid, behave as boid, else, random walk
        afraid = false 
    end
    
    methods
        
    end
    
end
