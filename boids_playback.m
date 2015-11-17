%% Playback Engine %%

% Playback engine to replay the exact scenario.
% Load a "history" variable by loading one of the trainset data: e.g. trainset1.mat
% Speed of playback can be controlled via the parameters.

%% Blank slate
clc
fprintf('Script executing...\n');
close all
clearvars -except history

%% Playback Parameters

fps = 30;
playback_speed = 2; % >0x to 4x (decimals accepted)


%% Simulation Parameters (please edit to match playback data)
field_size = 500; % edit if required

%% Initialisation of Playback

% initialise the figure
fig1 = figure(1);
hold on
axis([0 field_size 0 field_size])
plotHandle = plot(0,0,'ok');
sheepdogHandle = plot(0,0,'xr');

i = 1;
prev_toc = 0;
tic
while i <= size(history.sheep_x,1)
    if toc - prev_toc >= 1/fps
        plotHandle.set('XData', history.sheep_x(i,:), 'YData', history.sheep_y(i,:));
        sheepdogHandle.set('XData',history.mouse_pos(i,1),'YData',history.mouse_pos(i,2));
        prev_toc = toc;
    end
    i = i + 1;
    pause(0.04/playback_speed)
end