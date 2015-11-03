function com = boids_com( boids_array )
%FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here
    sum = [0,0];
    for i = 1:numel(boids_array)
        sum = sum + boids_array(i).position;
    end
    com = sum / numel(boids_array);
end





