%% Sheepdog Trainer

% Trains the sheepdog neural network using the datasets that you have

%% Parameters

% names of the files training sessions
filenames = ['trainset01';'trainset02';'trainset03';'trainset04';'trainset05'; ...
    'trainset06';'trainset07';'trainset08';'trainset09';'trainset10'];

learn_rate = 1e-6; % 5e-8 NN, 2e-7 WNN, 1e-6? DNN
epochs = 10;

%% Prepare Training Data

input_set = [];
output_set = [];

for i = 1:size(filenames,1)
    load(filenames(i,:));
    input_set = [input_set; history.sheep_x history.sheep_y history.mouse_pos];
    input_set(end,:) = [];
    output_set = [output_set; history.mouse_velocity];
end

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
    if mod(cycle,2) == 0 % plot every x number of epochs
        avg_error = TestNN(NN,input_set,output_set);
        fprintf('\n         Test error (training set) = %.6f\n\n',avg_error);
        plot(cycle,avg_error, '.k')
        drawnow
        if avg_error < 1e-4
            fprintf('\nFully converged!\n');
            break
        end
        if isnan(avg_error)
            error('Training diverged.');
        end
    end
end
hold off