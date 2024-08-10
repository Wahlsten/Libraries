clear
clc
error     = error_obj;
error.N   = [20 30 40 50; 20 30 40 50]; %[20, 30, 40, 50; 20, 30, 40, 50];
error.Acc = [2, 4, 6, 8];
problem   = 'Advection2D'; % 'Advection2D', 'Navier-Stokes'
data      = 'sin2pi'; % 'sincos_smooth', 'sin2pi', 'sin'
error     = initialize(error, problem, data);