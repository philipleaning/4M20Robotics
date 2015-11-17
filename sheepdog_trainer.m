%% Sheepdog Trainer

% Trains the sheepdog neural network using the datasets that you have

%% Parameters

% names of the files training sessions
filenames = ['trainset01';'trainset02';'trainset03';'trainset04';'trainset05'; ...
    'trainset06';'trainset07';'trainset08';'trainset09';'trainset10';'trainset11'];

learn_rate = 8e-4; % 1e-4
epochs = 30; % 1epoch = 10secs, 6ep/min 360ep/hr

%% Prepare Training Data

input_set = [];
output_set = [];

for i = 1:size(filenames,1)
    load(filenames(i,:));
    input_set = [input_set; history.grid_status history.mouse_pos./500]; % note the normalisation of the input
    input_set(end,:) = []; % take away last row since mouse velocity has one less row
    output_set = [output_set; history.mouse_velocity];
end

% normalisation to allow network to operate (this must be un-normalised
% during actual playback)
output_set = output_set + 50;
output_set = output_set./100;

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