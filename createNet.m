
%% Load Training Data
% Load '.mat' data is located in
fileNames = ['TrainingData'];
load(fileNames);
% Store data in local variables in correct format (each step as collumns)
input = [];
output = [];
input = [input inputDataForNet'];
output = [output outputDataForNet'];

%% Create and Train Net
% Create net
net = feedforwardnet(10);
% Train net using LM
[net,tr] = trainlm(net,input,output);

%% Test & Save Net
% Test
y = net(inputDataForNet');
perf = perform(net,y,outputDataForNet')
% Save in '.mat' file
save('TrainedNet', 'net');