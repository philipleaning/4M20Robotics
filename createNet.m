
fileNames = ['TrainingData'];

load(fileNames);

input = [];
output = [];
input = [input inputDataForNet'];
output = [output outputDataForNet'];

net = feedforwardnet(10);
[net,tr] = trainlm(net,input,output);

y = net(inputDataForNet');

perf = perform(net,y,outputDataForNet')

save('TrainedNet', 'net');