
%% Load Training Data
% Load '.mat' data is located in
fileNames = ['TrainingData'];
load(fileNames);
% Store data in local variables in correct format (each step as collumns)
input = [];
output = [];
input = [input inputDataForNet'];
output = [output outputDataForNet'];

x = input
t = output;

%% Net Creation and Training
% Choose a Training Function
% For a list of all training functions type: help nntrain
% 'trainlm' is usually fastest.
% 'trainbr' takes longer but may be better for challenging problems.
% 'trainscg' uses less memory. Suitable in low memory situations.
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

% Create a Fitting Network
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize,trainFcn);

% Setup Division of Data for Training, Validation, Testing
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
y = net(x);
e = gsubtract(t,y);
performance = perform(net,t,y)

% View the Network
%view(net);

% Save in '.mat' file
save('TrainedNet', 'net');
