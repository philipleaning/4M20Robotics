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
number_of_boids = 10;
boids_array = Boid.empty;
max_speed = 5;

%% Initialisation
fprintf('Initialising simulation...');

% initialise the figure
fig1 = figure(1);
hold on
axis([0 field_size 0 field_size])
plotHandle = plot(0,0, 'o');
dogPlot = plot(0,0,'x');
dog_1 = Dog();

% Load neural net
%load('TrainedNet', 'net');
load TrainingData5RoundsOf200


% Create nets

%nets = net.empty()

net1 = createNetFromFileName('TrainingData001GPS_50RoundsOf200');
net2 = createNetFromFileName('TrainingData001GPS_30RoundsOf200');
net3 = createNetFromFileName('TrainingData001GPS_20RoundsOf200');
net4 = createNetFromFileName('TrainingData001GPS_10RoundsOf200');
net5 = createNetFromFileName('TrainingData001GPS_5RoundsOf200');
net6 = createNetFromFileName('TrainingData001GPS_2RoundsOf200');
net7 = createNetFromFileName('TrainingData001GPS_1RoundsOf200');

% Stick in array
nets = {net1 net2 net3 net4 net5 net6 net7};

% Store end variance of each net
endVarianceForEachNet = [];

% create boids
for i = 1:number_of_boids
    boids_array(end+1) = Boid();
    boids_array(i).position = [rand*field_size, rand*field_size];
end

%% Simulation
fprintf('Running simulation...')


for n=1:numel(nets)
    n
    net = nets{n};

    % Store variance of each trial
    trialVariances = [];
    
    for i=1:1000
        if mod(i,100) == 0 
            i
        end
        %% Reload shob and positions every 500 steps. Store variance of shobs
        if mod(i,200) == 0 
            total_pos = 0;
            for x = 1:number_of_boids
                total_pos = total_pos + boids_array(x).position;
            end
            cenOfMass = total_pos/number_of_boids;
            sumVar = 0;
            for y = 1:number_of_boids
                a=(boids_array(y).position-cenOfMass);
                boidDist=sqrt(a(1)*a(1)+a(2)*a(2));
                var=boidDist*boidDist;
                sumVar = sumVar + var;
            end
            currentVariance = sumVar/number_of_boids;
            for k = 1:number_of_boids
                boids_array(k).position = [rand*field_size, rand*field_size];
                boids_array(k).velocity = [0,0];
            end 
            dog_1.position = [rand*field_size, rand*field_size];
            dog_1.velocity = [0,0];
            trialVariances = [trialVariances; currentVariance];
        end

        %% Count shobs in quadrants, update velocity with net
        % For each boid, see if close enough to sheepdog (mouse) to be afraid
        for j = 1:numel(boids_array)
            boid = boids_array(j);      
            distance_from_dog = norm(boid.position - dog_1.position);
            if distance_from_dog < 100 
               boid.afraid = true;
            else
               boid.afraid = false;
            end
        end

        % Reset sheep mass count to 0
        dog_1.sheepMass = [0 0 0 0];
        sheepPos = [];
        % Count sheep in quadrants
        for x = 1:numel(boids_array)
           b = boids_array(x);
           sheepPos = [sheepPos;b.position];
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

        %% Apply Boid rules to scared Shobs
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
                    v3 = 0;
                end
                % FEAR!!: If afraid add velocity away from dog
                if boid.afraid 
                   fear_velocity = (boid.position - dog_1.position)/norm(boid.position - dog_1.position) * 0.8; 
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
        
        %% Bound boids approximately to field, bounce them back with a velocity
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
        
        %% Move Dog and Bound Dog to field
        if i > 20
            %input = [dog_1.sheepMass]; %for quad
            input = [sheepPos(:,1);sheepPos(:,2);dog_1.position']; %for GPS
            %dog_1.velocity = net(input')'; %for quad
            dog_1.velocity = net(input)'; %for GPS
        end
        dog_1.sheepMassHistory = [dog_1.sheepMassHistory; dog_1.sheepMass];

        speed = norm(dog_1.velocity);
        if speed > 16
              dog_1.velocity = (dog_1.velocity / speed) * 16;
        end
        dog_1.velocity;
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
        
        %% Update figures
        %plotHandle.set('XData', boids_x_pos, 'YData', boids_y_pos);
        %dogPlot.set('XData', dog_1.position(1), 'YData', dog_1.position(2));
        %pause(0.05);
    end
    %% Average End Variance and Store
    n
    averageVariance = mean(trialVariances);
    endVarianceForEachNet = [endVarianceForEachNet; averageVariance]
    
end
trainLength = [50 30 20 10 5 2 1];
plot(trainLength,endVarianceForEachNet)
