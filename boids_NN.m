%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%          Sheep Simulator                                      %%%%%%
%%%%%%          4M20 Robotics, coursework template                   %%%%%%
%%%%%%          University of Cambridge                              %%%%%%
%%%%%%          Michaelmas 2015                                      %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sheep Simulation                                                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Blank slate
clc
fprintf('Script executing...\n');
close all
clearvars -except NN

%% Parameters
field_size = 500;
number_of_boids = 10;
boids_array = Boid.empty;
max_speed = 5;  % boid maximum speed

% Bottom left cage coordinates  which stops simulation once all boids are
% inside. Top right coordinates are 500,500
cage_x = 450; 
cage_y = 450;

%% Initialisation

fprintf('Initialising simulation...\n');

% initialise the figure
fig1 = figure(1);
hold on
axis([0 field_size 0 field_size])
plotHandle = plot(0,0,'ok');
sheepdogHandle = plot(0,0,'xr');

% create boids
boids_x_pos = zeros(1,number_of_boids);
boids_y_pos = zeros(1,number_of_boids);

for i = 1:number_of_boids
    boids_array(end+1) = Boid();
    boids_array(i).position = [rand*field_size, rand*field_size];
    boids_x_pos(i) = boids_array(i).position(1);
    boids_y_pos(i) = boids_array(i).position(2);
end

% initial starting point of sheepdog 
mousePoint = [50 50];
%% Simulation
fprintf('Running simulation...\n')

for i=1:10000
    
    NN = RunNN(NN,[boids_x_pos boids_y_pos mousePoint]);
    NN.output
    mousePoint = mousePoint + NN.output; % Get sheepdog position.
    
    % stop the sheepdog from exiting the field
    mousePoint(mousePoint > 500) = 500;
    mousePoint(mousePoint < 0) = 0;
    
    % Handle each boid at a time
    for j = 1:numel(boids_array)
        
        boid = boids_array(j); % The current boid
        
        distance_from_dog = norm(boid.position - mousePoint);
        
        % For each boid, see if close enough to sheepdog (mouse) to be afraid
        
        if distance_from_dog < 100 
            boid.afraid = true;

            % Use boid rules when afraid %
            
            % For each boid: sum the vectors from applying the 3 rules
            
            % Boids move to centre of gravity, 1% a step
            pos_com = boids_com(boids_array);
            v1 = ( pos_com - boid.position ) / 300;

            % Boids avoid each other, repel from all boids within 3 units
            v2 = [0,0];
            for k = 1:numel(boids_array)
               if boids_array(k) ~= boid % if not current boid
                   if norm(boids_array(k).position - boid.position) < 3
                        v2 = v2 - (boids_array(k).position - boid.position);
                   end
               end
            end

            % Boids tend to each other's velocities
            % Add 1/8th of difference between current boid velocity and perceived velocity of the swarm
            v3 = (boid_perceived_vel(boids_array, boids_array(j)) - boid.velocity) / 8;

            % FEAR!!: If afraid add velocity away from dog
            fear_velocity = (boid.position - mousePoint)/norm(boid.position - mousePoint) * 0.8;

            % Update boid velocity by summing contribution from each rule
            boid.velocity = boid.velocity + v1 + v2 + v3 + fear_velocity;

            % Limit velocity
            speed = norm(boid.velocity);
            if speed > max_speed
                boid.velocity = (boid.velocity / speed) * max_speed;
            end

            % Update boids position
            boid.position = boid.position + boid.velocity; 

        else
            boid.afraid = false;

            % Random walk %
            
            % Update boids velocity, random walk
            delta_velocity = [(randn)*0.1,(randn)*0.1];
            boid.velocity = 0.93*boid.velocity + delta_velocity;
            boid.position = boid.position + boid.velocity;
        end
        
        % Bound boids approximately to field, bounce them back with a velocity change if they leave
        if boid.position(1) < 30 && boid.velocity(1) < 0
           boid.velocity(1) = boid.velocity(1) + 2;
        end
        if field_size - boid.position(1) < 30 && boid.velocity(1) > 0
           boid.velocity(1) = boid.velocity(1) - 2;
        end
        if boid.position(2) < 30 && boid.velocity(2) < 0
           boid.velocity(2) = boid.velocity(2) + 2;
        end
        if field_size - boid.position(2) < 30 && boid.velocity(2) > 0
           boid.velocity(2) = boid.velocity(2) - 2;
        end
        
        % Update boids position arrays
        boids_x_pos(j) = boid.position(1);
        boids_y_pos(j) = boid.position(2);
        
    end
    
    plotHandle.set('XData', boids_x_pos, 'YData', boids_y_pos)
    sheepdogHandle.set('XData',mousePoint(1),'YData',mousePoint(2))
    
    % Stop simulation once boids have reached above x,y because
    % they have been successfully herded
    if sum(boids_x_pos <= cage_x) == 0
        if sum(boids_y_pos <= cage_y) == 0
            break
        end
    end
    
    pause(0.01);
end