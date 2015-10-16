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
clear all
close all
clc

%% Parameters

no_of_sheep = 50;
dt = 1e-2;          % simulation step size (smallest is 1e-2)
sim_length = 30;    % simulation length (no. of seconds)
sheep_mvt = 0.1;      % amount of sheep movement per simulation step
fps = 25;           % update rate
start_mvt = 0.008;   % chance of sheep to start movement
cont_mvt = 0.95;     % chance of sheep to continue moving once moving
deg_of_freedom = 30;% the arc which the sheep is free to choose to head in

%% Initialisation

% initialise some variables
sheep_pos_x = zeros(sim_length/dt,no_of_sheep); % preallocate
sheep_pos_y = zeros(sim_length/dt,no_of_sheep); % preallocate
sheep_move = zeros(1,no_of_sheep);              % preallocate
prev_angle = zeros(1,no_of_sheep);              % preallocate
prev_toc = 0;

% initialise the figure
figure(1)
hold on
axis([0 100 0 100])

for i = 1:no_of_sheep
    sheep_pos_x(1,i) = rand*60 + 20; % randomise from 20-80
    sheep_pos_y(1,i) = rand*60 + 20; % randomise from 20-80
    prev_angle(1,i) = rand*2*pi;
    sheep(i) = line(sheep_pos_x(1,i),sheep_pos_y(1,i),'color',[0 0 0],'Marker','.','MarkerSize',10); 
end

%% Simulation

tic
for i = 2:sim_length/dt
    for sheep_no = 1:no_of_sheep
        if sheep_move(1,sheep_no) == 1
            if rand < cont_mvt
                % move the sheep
                angle = (rand-0.5)*deg_of_freedom/180*pi + prev_angle(1,sheep_no);
                prev_angle(1,sheep_no) = angle;
                sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no) + cos(angle)*sheep_mvt;
                sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no) + sin(angle)*sheep_mvt;
                % prevent sheep from exceeding bounds
                if abs(sheep_pos_x(i,sheep_no)-50) > 50
                   sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no);
                   prev_angle(1,sheep_no) = rand*2*pi; % also reset the angle to random one so the sheep doesn't get stuck
                end
                if abs(sheep_pos_y(i,sheep_no)-50) > 50
                   sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no);
                   prev_angle(1,sheep_no) = rand*2*pi; % also reset the angle to random one so the sheep doesn't get stuck
                end
            else
                sheep_move(1,sheep_no) = 0;
                % don't move sheep
                sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no);
                sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no);
            end
        else
            if rand < start_mvt
                sheep_move(1,sheep_no) = 1;
                % move the sheep
                angle = (rand-0.5)*deg_of_freedom/180*pi + prev_angle(1,sheep_no);
                prev_angle(1,sheep_no) = angle;
                sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no) + cos(angle)*sheep_mvt;
                sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no) + sin(angle)*sheep_mvt;
                % prevent sheep from exceeding bounds
                if abs(sheep_pos_x(i,sheep_no)-50) > 50
                   sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no);
                   prev_angle(1,sheep_no) = rand*2*pi; % also reset the angle to random one so the sheep doesn't get stuck
                end
                if abs(sheep_pos_y(i,sheep_no)-50) > 50
                   sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no); 
                   prev_angle(1,sheep_no) = rand*2*pi; % also reset the angle to random one so the sheep doesn't get stuck
                end
            else
                % don't move sheep
                sheep_pos_x(i,sheep_no) = sheep_pos_x(i-1,sheep_no);
                sheep_pos_y(i,sheep_no) = sheep_pos_y(i-1,sheep_no);
            end
        end
    end
    now_toc = toc;
    if now_toc - prev_toc >= 1/fps
        for j = 1:no_of_sheep
            set(sheep(j),'xdata',sheep_pos_x(i,j),'ydata',sheep_pos_y(i,j));
        end
        drawnow
        prev_toc = now_toc;
    end
    pause(dt)
end