function perceived_velocity = boid_perceived_vel( boids_array, current_boid )
% BOID_VELOCITY 
% Finds average velocity of all boids except current, the "perceived velocity"
%
    perceived_velocity = [0,0];

    for i=1:numel(boids_array) 
       if boids_array(i) ~= current_boid
          perceived_velocity = perceived_velocity + boids_array(i).velocity;
       end
    end
    
    perceived_velocity = perceived_velocity / (numel(boids_array) - 1);
    
end

