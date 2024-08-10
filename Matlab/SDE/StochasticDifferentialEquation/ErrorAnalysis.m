clear
clc

%% Strong Convergence

errorS         = error_obj;
errorS.N       = 2^9;
errorS.M       = 1000;
errorS.R       = 2.^(0:4);
errorS.method  = {'Euler', 'Milstein'};
errorS.NoExact = 0; % No exact solution: no = 0, yes = 1
errorS.Acc     = [1/2, 1];

problemS       = 'Black-Scholes'; % 'Black-Scholes2'
dataS          = 'constant';
errorS         = initialize(errorS, problemS, dataS);
errorS         = StrongConvergence(errorS, problemS, dataS);

%% Weak Convergence

randn('state', 100)
errorW         = error_obj;
errorW.N       = 2.^(5:9);
errorW.M       = 50000;
errorW.R       = 1;
errorW.method  = {'Euler', 'Milstein'}; % 'Euler', 'Milstein', 'RK'
errorW.NoExact = 0; % No exact solution: no = 0, yes = 1
errorW.Acc     = [1, 1];

problemW       = 'Black-Scholes2'; % 'Black-Scholes2'
dataW          = 'constant';
errorW         = initialize(errorW, problemW, dataW);
errorW         = WeakConvergence(errorW, problemW, dataW);