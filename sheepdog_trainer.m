%% Sheepdog Trainer

% Trains the sheepdog neural network using the datasets that you have

%% Parameters

% Type of training: training or reinforcement
training_type = 'training';

learn_rate = 2e-3; % 5e-3
epochs = floor(120/1.2); % 1epoch = ~50secs, 1.2ep/min 72ep/hr

%% Prepare Training Data

% prepare the file names for reading
if strcmp(training_type,'training')
    stem = 'trainset';
    filenames = ['trainset01.mat'];
elseif strcmp(training_type,'reinforcement')
    stem = 'success';
    filenames = ['success01.mat'];
else
    error('no/wrong training scheme selected.')
end

i = 2;
while 1
    if i < 10
        filenames(i,:) = strcat(stem,'0',num2str(i),'.mat');
    else
        filenames(i,:) = strcat(stem,num2str(i),'.mat');
    end
    if exist(filenames(i,:),'file')
        i = i + 1;
    else
        filenames(i,:) = [];
        break
    end
end



input_set = [];
output_set = [];

for i = 1:size(filenames,1)
    load(filenames(i,:));
    input_set = [input_set; history.sheep_x history.sheep_y history.mouse_pos];
    input_set(end,:) = []; % take away last row since mouse velocity has one less row
    output_set = [output_set; history.mouse_velocity];
end

% normalisation to allow network to operate (this must be un-normalised
% during actual playback)
output_set = output_set + 50;
output_set = output_set./100;
input_set = input_set./500;

clear history

%% Prepare to train the network %%

close all
figure(1)
xlim([0 epochs])
hold on

% Resetting the card prevents some weird slowdown problems
g_card = gpuDevice(1);
reset(g_card);

%% Training Routine %%

for cycle = 1:epochs
    NN = TrainNNGPU(NN,input_set,output_set,learn_rate);
    fprintf('Epochs completed: %d\n',cycle);
    if mod(cycle,1) == 0 % plot every x number of epochs
        avg_error = TestNN(NN,input_set,output_set);
        fprintf('\n         Avg test error^2 (training set) = %.6f\n\n',avg_error);
        plot(cycle,avg_error, '.k')
        drawnow
        if avg_error < 1e-6
            fprintf('\nFully converged!\n');
            break
        end
        if isnan(avg_error)
            error('Training diverged.');
        end
    end
end
hold off

avg_error = TestNN(NN,input_set,output_set);
fprintf('\n***** Actual average error = %.6f *****\n\n',sqrt(avg_error)*100);