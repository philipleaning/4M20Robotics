
fileNames = ['TrainingData'];

load(fileNames);

input = [];
output = [];
input = [input inputDataForNet'];
output = [output outputDataForNet'];

net = fitnet(10);
[net,tr] = train(net,input,output);

y = net(inputDataForNet');

perf = perform(net,y,outputDataForNet')