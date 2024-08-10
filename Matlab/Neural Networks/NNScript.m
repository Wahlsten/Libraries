%% Data
[x,t] = simplenarx_dataset;


%% Parameters
sizeNet = [25, 25];
numInputDelays = 10;
numOutputDelays = 10;

%% Network
net = narxnet(1:numInputDelays, 1:numOutputDelays, sizeNet);
[X,Xi,Ai,T] = preparets(net,x,{},t);
net = train(net,X,T,Xi,Ai);
view(net)
Y = net(X,Xi,Ai);
perf = perform(net,Y,T);