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

% Simulation Options
sim_length = 2000;    % simulation length (100 is very roughly 1sec of playback at 1x speed)
no_of_sheep = 50;     % self explanatory (cannot have too many or placement will fail!)
field_size = 100;     % size of the field that the sheep are grazing on
sheep_mvt = 0.5;      % amount of sheep movement per simulation step
p_start_mvt = 0.008;  % chance of sheep to start movement
p_cont_mvt = 0.90;    % chance of sheep to continue moving once moving
deg_of_freedom = 30;  % the arc which the sheep is free to choose to head in
sheep_dist = 4;       % min distance between adjacent sheep


% Playback Options
fps = 25;             % draw rate (usually around 25 gives a smooth playback)
plot_paths = false;   % choose to plot paths of sheep? (very resource intensive! reduce draw rate, e.g. 5fps to help!)
playback_speed = 1;   % playback speed for animation. 1x, 2x, 5x, 10x, etc.


%% Initialisation
fprintf('Initialising simulation...');

% initialise some variables
sheep_pos_x = zeros(sim_length,no_of_sheep); % position x of sheep
sheep_pos_y = zeros(sim_length,no_of_sheep); % position y of sheep
sheep_move = zeros(1,no_of_sheep);           % elements are true for sheep which are moving
prev_angle = zeros(1,no_of_sheep);           % previous heading of sheep
prev_toc = 0;
prev_i = 1;

% initialise the figure
figure(1)
hold on
axis([0 field_size 0 field_size])

% initialise positions of sheep. Algorithm also watches out and doesn't
% place sheep too close to each other.
tic
for i = 1:no_of_sheep
    sheep_pos_x(1,i) = rand*(field_size*0.8) + field_size*0.1; % initialise at middle 80% of the field
    sheep_pos_y(1,i) = rand*(field_size*0.8) + field_size*0.1; % initialise at middle 80% of the field
    k = 1;
    % check distance between current sheep being placed and other sheep
    % before it
    while k <= i
        if k ~= i
            if sqrt(abs(sheep_pos_x(1,i)-sheep_pos_x(1,k))^2 + abs(sheep_pos_y(1,i)-sheep_pos_y(1,k))^2) < sheep_dist
                while sqrt(abs(sheep_pos_x(1,i)-sheep_pos_x(1,k))^2 + abs(sheep_pos_y(1,i)-sheep_pos_y(1,k))^2) < sheep_dist
                    sheep_pos_x(1,i) = rand*(field_size*0.8) + field_size*0.1; % initialise at middle 80% of the field
                    sheep_pos_y(1,i) = rand*(field_size*0.8) + field_size*0.1; % initialise at middle 80% of the field
                    if toc > 2
                        error('Cannot find initial placement for 1 or more sheep. Likely too many sheep!')
                    end
                end
                k = 0;
            end
        end
        k = k + 1;
    end
    prev_angle(1,i) = rand*2*pi; % randomise the direction sheep face
    sheep(i) = line(sheep_pos_x(1,i),sheep_pos_y(1,i),'color',[0 0 0],'Marker','.','MarkerSize',15); % marker properties to visualise sheep 
end


fprintf('done\n')
%% Simulation
fprintf('Running simulation...')

tic
for i = 2:sim_length
    for sheep_no = 1:no_of_sheep
        if sheep_move(1,sheep_no) == true
            if rand < p_cont_mvt
                % do nothing (code to implement movement executed later)
            else
                sheep_move(1,sheep_no) = false;
            end
        else
            if rand < p_start_mvt
                sheep_move(1,sheep_no) = true;
            else
                % do nothing (code to implement non-movement executed later)
            end
        end
        % code to move sheep
        if sheep_move(1,sheep_no) == true
            angle = (rand-0.5)*deg_of_freedom/180*pi + prev_angle(1,sheep_no);
            prev_angle(1,sheep_no) = angle;
            sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no) + cos(angle)*sheep_mvt;
            sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no) + sin(angle)*sheep_mvt;
            
            % prevent sheep from walking into each other
            k = 1;
            while k <= no_of_sheep
                if k ~= sheep_no
                    if sqrt(abs(sheep_pos_x(i,sheep_no)-sheep_pos_x(i-1,k))^2 + abs(sheep_pos_y(i,sheep_no)-sheep_pos_y(i-1,k))^2) < sheep_dist
                        while sqrt(abs(sheep_pos_x(i,sheep_no)-sheep_pos_x(i-1,k))^2 + abs(sheep_pos_y(i,sheep_no)-sheep_pos_y(i-1,k))^2) < sheep_dist
                            angle = rand*2*pi; % randomise the new angle it should face and attempt to walk
                            prev_angle(1,sheep_no) = angle;
                            sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no) + cos(angle)*sheep_mvt*1.1;
                            sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no) + sin(angle)*sheep_mvt*1.1;
                        end
                        k = 0;
                    end
                end
                k = k + 1;
            end
            % prevent sheep from exceeding bounds
            if abs(sheep_pos_x(i,sheep_no)-field_size/2) > field_size/2
               sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no);
               prev_angle(1,sheep_no) = rand*2*pi; % also reset the angle to random one so the sheep doesn't get stuck
            end
            if abs(sheep_pos_y(i,sheep_no)-field_size/2) > field_size/2
               sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no);
               prev_angle(1,sheep_no) = rand*2*pi; % also reset the angle to random one so the sheep doesn't get stuck
            end
        
        % code to not move sheep
        else
            sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no);
            sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no);
        end
    end
    % fprintf('Step: %d\n',i);
end

fprintf('done\n');
fprintf('Simulation took %.3f seconds.\n',toc);
%% Playback
fprintf('Begin playback...');
tic
for i = 2:playback_speed:sim_length
    now_toc = toc;
    if now_toc - prev_toc >= 1/fps
        for j = 1:no_of_sheep
            set(sheep(j),'xdata',sheep_pos_x(i,j),'ydata',sheep_pos_y(i,j));
            if plot_paths
                if (sheep_pos_x(prev_i,j) ~= sheep_pos_x(i,j)) || (sheep_pos_y(prev_i,j) ~= sheep_pos_y(i,j))
                    plot([sheep_pos_x(prev_i,j),sheep_pos_x(i,j)],[sheep_pos_y(prev_i,j),sheep_pos_y(i,j)]);
                end
            end
        end
        drawnow
        prev_i = i;
        prev_toc = now_toc;
    end
    pause(0.001)
end
fprintf('done\n');
fprintf('Script finished executing.\n');
