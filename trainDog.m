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
field_size = 500;
number_of_boids = 15;
boids_array = Boid.empty;
max_speed = 5;

%% Stored data for training NN


% Playback Options
fps = 25;             % draw rate (usually around 25 gives a smooth playback)
plot_paths = false;   % choose to plot paths of sheep? (very resource intensive! reduce draw rate, e.g. 5fps to help!)
playback_speed = 1;   % playback speed for animation. 1x, 2x, 5x, 10x, etc.


%% Initialisation

fprintf('Initialising simulation...');

set (gcf, 'WindowButtonMotionFcn', @mouseMove);

% initialise the figure
fig1 = figure(1);
hold on
axis([0 field_size 0 field_size])
plotHandle = plot(0,0, 'o');
dogPlot = plot(0,0,'x');
dog_1 = Dog();
inputData = [0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0 ];
inputData2 = zeros(5,5)
inputData3 = zeros(5,20)
sheepMassHistory2 = zeros(2000,4);
% create boids
for i = 1:number_of_boids
    boids_array(end+1) = Boid();
    boids_array(i).position = [rand*field_size, rand*field_size];
end

%% Simulation
fprintf('Running simulation...')

for i=1:100
    
    if mod(i,200)==0 
        for k = 1:number_of_boids
            boids_array(k).position = [rand*field_size, rand*field_size];
        end 
        
    end
    
    % For each boid, see if close enough to sheepdog (mouse) to be afraid
    for j = 1:numel(boids_array)
        boid = boids_array(j);
        % Get sheepdog (mouse pointer) position.
        pointMatrix = getMousePoint;
      
        if ~isempty(pointMatrix)
            mousePoint = pointMatrix(1,1:2);
            dog_1.position = mousePoint;     
            distance_from_dog = norm(boid.position - dog_1.position);
            if distance_from_dog < 100 
               boid.afraid = true;
            else
               boid.afraid = false;
            end
        end 
        
    end
    if ~isempty(pointMatrix)
        pointMatrix = getMousePoint;
        mousePoint = pointMatrix(1,1:2);
        
        dog_1.position = mousePoint;     
        dog_1.velocity = dog_1.position - dog_1.positionHistory(max(i-1,1),:);        
        
        dog_1.velocityHistory = [dog_1.velocityHistory; dog_1.velocity];
        dog_1.positionHistory = [dog_1.positionHistory; dog_1.position];
        
        % Reset sheep mass count to 0
        dog_1.sheepMass = [0 0 0 0];
        % Count sheep in quadrants
        for x = 1:numel(boids_array)
           dog_1.sheepPosHistory(i,2*x-1)= dog_1.position(1)-b.position(1); 
           dog_1.sheepPosHistory(i,2*x)= dog_1.position(2)-b.position(2);
           b = boids_array(x);
           distMass = (1/norm(dog_1.position-b.position));
      
           if b.position(1) < dog_1.position(1) && b.position(2) < dog_1.position(2) %bot left
               
               dog_1.sheepMass(1)=dog_1.sheepMass(1)+distMass;
           end
           if b.position(1) > dog_1.position(1) && b.position(2) < dog_1.position(2) % bot right
               dog_1.sheepMass(2)=dog_1.sheepMass(2)+distMass;
           end
           if b.position(1) < dog_1.position(1) && b.position(2) > dog_1.position(2) % top left
               dog_1.sheepMass(3)=dog_1.sheepMass(3)+distMass;
           end
           if b.position(1) > dog_1.position(1) && b.position(2) > dog_1.position(2) % top right
               dog_1.sheepMass(4)=dog_1.sheepMass(4)+distMass;
           end
        end
        dog_1.sheepMassHistory = [dog_1.sheepMassHistory; dog_1.sheepMass];
    end
    
     %{
    t = cputime;
    surf(peaks(40));
    e = cputime-t
    %}
    % For each boid: sum the vectors from applying the 3 rules
    for j = 1:numel(boids_array)
        % The current boid
        boid = boids_array(j);
        
        % If afraid use boid rules, else, random walk
        if boid.afraid == true
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

            % Boids tend to each other's velocities, 
            % Add 1/8th of difference between current boid velocity and perceived velocity of the swarm
            if number_of_boids > 1
                v3 = (boid_perceived_vel(boids_array, boids_array(j)) - boid.velocity) / 8;
            else
                v3 = 0
            end
            % FEAR!!: If afraid add velocity away from dog
            if boid.afraid 
               fear_velocity = (boid.position - mousePoint)/norm(boid.position - mousePoint) * 0.8; 
            end
            % Update boid velocity by summing contribution from each rule
            boid.velocity = boid.velocity + v1 + v2 + v3 + fear_velocity;
            % Limit velocity
            speed = norm(boid.velocity);
            if speed > max_speed
                boid.velocity = (boid.velocity / speed) * max_speed;
            end
            % Update boids position
            boid.position = boid.position + boid.velocity; 
        else % Random walk
            % Update boids velocity, random walk
            delta_velocity = [(randn)*0.1,(randn)*0.1];
            boid.velocity = 0.93*boid.velocity + delta_velocity;
            boid.position = boid.position + boid.velocity;
        end
       
    end

    % Bound boids approximately to field, bounce them back with a velocity
    % change if they leave
    for j=1:numel(boids_array)
       b = boids_array(j);
       if b.position(1) < 30 && b.velocity(1) < 0
           b.velocity(1) = b.velocity(1) + 2;
       end
       if field_size - b.position(1) < 30 && b.velocity(1) > 0
           b.velocity(1) = b.velocity(1) - 2;
       end
       if b.position(2) < 30 && b.velocity(2) < 0
           b.velocity(2) = b.velocity(2) + 2;
       end
       if field_size - b.position(2) < 30 && b.velocity(2) > 0
           b.velocity(2) = b.velocity(2) - 2;
       end
    end
    boids_x_pos = zeros(1,number_of_boids);
    boids_y_pos = zeros(1,number_of_boids);
    
    for j = 1:numel(boids_array)
        b = boids_array(j); % Current boid
        boids_x_pos(j) = b.position(1);
        boids_y_pos(j) = b.position(2);
    end
      
    acc = norm(dog_1.deltaVelocity);
    if acc > 4
          dog_1.deltaVelocity = (dog_1.deltaVelocity / acc) * 7;
    end
    
    speed = norm(dog_1.velocity);
    if speed > 200
          dog_1.velocity = (dog_1.velocity / speed) * 16;
    end
    dog_1.position = dog_1.position + dog_1.velocity;
    
    if dog_1.position(1) < 30 && dog_1.velocity(1) < 0
        dog_1.position(1) = 30;
          %dog_1.velocity(1) = dog_1.velocity(1) + 15;
    end
    if field_size - dog_1.position(1) < 30 && dog_1.velocity(1) > 0
        dog_1.position(1) = field_size - 30;
          %dog_1.velocity(1) = dog_1.velocity(1) - 15;
    end
    if dog_1.position(2) < 30 && dog_1.velocity(2) < 0   
        dog_1.position(2) = 30;
          %dog_1.velocity(2) = dog_1.velocity(2) + 15;     
    end       
    if field_size - dog_1.position(2) < 30 && dog_1.velocity(2) > 0
        dog_1.position(2) = field_size - 30;
          %dog_1.velocity(2) = dog_1.velocity(2) - 15;      
    end                       
    
    
    plotHandle.set('XData', boids_x_pos, 'YData', boids_y_pos)
    dogPlot.set('XData', dog_1.position(1), 'YData', dog_1.position(2))
    pause(0.05);
end

%% Save Data

inputDataForNet = dog_1.sheepMassHistory;
outputDataForNet = dog_1.velocityHistory;

save('TrainingData', 'inputDataForNet', 'outputDataForNet')