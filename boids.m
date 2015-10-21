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
clear all
close all

%% Parameters
field_size = 300;
number_of_boids = 20;
boids_array = Boid.empty;

% Playback Options
fps = 25;             % draw rate (usually around 25 gives a smooth playback)
plot_paths = false;   % choose to plot paths of sheep? (very resource intensive! reduce draw rate, e.g. 5fps to help!)
playback_speed = 1;   % playback speed for animation. 1x, 2x, 5x, 10x, etc.


%% Initialisation

fprintf('Initialising simulation...');

% initialise the figure
fig1 = figure(1);
hold on
axis([0 field_size 0 field_size])
plotHandle = plot(0,0, 'o');

% create boids
for i = 1:number_of_boids
    boids_array(end+1) = Boid();
    boids_array(i).position = [rand*100, rand*100];
end

%% Simulation
fprintf('Running simulation...')

for i=1:1000
    % For each boid: sum the vectors from applying the 3 rules
    for i = 1:numel(boids_array)
        % The current boid
        boid = boids_array(i);

        % Boids move to centre of gravity, 1% a step
        pos_com = boids_com(boids_array);
        v1 = ( pos_com - boid.position ) / 100;

        % Boids avoid each other, repel from all boids within 3 units
        v2 = [0,0];
        for j = 1:numel(boids_array)
           if boids_array(j) ~= boid % if not current boid
               if norm(boids_array(j).position - boid.position) < 3
                    v2 = v2 - (boids_array(j).position - boid.position);
               end
           end
        end

        % Boids tend to each other's velocities, 
        % Add 1/8th of difference between current boid velocity and perceived velocity of the swarm
        v3 = (boid_perceived_vel(boids_array, boids_array(i)) - boid.velocity) / 8;

        % Update boid velocity by summing contribution from each rule
        boid.velocity = boid.velocity + v1 + v2 + v3;
        % Update boids position
        boid.position = boid.position + boid.velocity;

    end

    % Bound boids approximately to field, bounce them back with a velocity
    % change if they leave
    for i=1:numel(boids_array)
       b = boids_array(i);
       if b.position(1) < 0 
           b.velocity(1) = 1;
       end
       if b.position(1) > field_size
           b.velocity(1) = -1;
       end
       if b.position(2) < 0
           b.velocity(2) = 1;
       end
       if b.position(2) > field_size
           b.velocity(2) = -1;
       end
    end
    boids_x_pos = zeros(1,number_of_boids);
    boids_y_pos = zeros(1,number_of_boids);
    
    for i = 1:numel(boids_array)
        b = boids_array(i); % Current boid
        boids_x_pos(i) = b.position(1);
        boids_y_pos(i) = b.position(2);
    end
    
    plotHandle.set('XData', boids_x_pos, 'YData', boids_y_pos)
    
    pause(0.05);
end