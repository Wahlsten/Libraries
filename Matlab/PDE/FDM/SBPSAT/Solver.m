function Solver(problemType, problem, data, acc, N, minv, maxv, geometry, plotSol)
% problemType: Problem type:       'AdvectionScalar1D', ..., 'HyperbolicSystem2D'
% problem:     Problem parameters: 'Advection1D', ..., 'Maxwell2D'
% data:        Data:               'sin', ..., 'constant'
% acc:         Accuracy:           2, 4, 6, 8
% N:           Number of points:   [Nx, ..., Nt]
% minv:        Lower limit:        [minx, ..., mint]
% maxv:        Upper limit:        [maxx, ..., maxt]
% geometry:    Geometry:           'Unit', 'Rectangle'
% plotSol:     Plot solution:      1 = 'yes', 0 = 'no'

eval(['addpath(''/Users/markuswahlsten/Desktop/Coding/Problem/', problemType, ''')'])

eval(['[sol, t] = ', problemType, 'Solver(', '''', problem, '''',',', num2str(acc),',', ...
    '''', data, '''', ',', '[', num2str(N), ']', ',', '[', num2str(minv) ']', ',', '[', ...
    num2str(maxv), ']', ',', '''', geometry, '''',');'])

eval(['rmpath(''/Users/markuswahlsten/Desktop/Coding/Problem/', problemType, ''')'])

fprintf(['======= Solving ', problemType, ' =======\n'])
fprintf('error = %d \n', sol.Solution.e_L2)
fprintf('time = %d \n', t)
    
if plotSol
    
    figure
    
    if length(N) == 2
        
        plot(sol.Solution.u_n(end, :))
        
    elseif length(N) == 3
    
        plot(sol.Solution.u_n(end, :))
        
    elseif length(N) == 4
        
    x = linspace(minv(1), maxv(1), N(1)+1);
    y = linspace(minv(2), maxv(2), N(2)+1);

    [X, Y] = meshgrid(x, y);

    for k = 1:size(sol.Solution.u_n, 2)

       u = reshape(sol.Solution.u_n(end - (N(1)+1) * (N(2)+1) + 1:end, k), N(1)+1, N(2)+1);
       surf(X,Y,u)
       axis([0 1 0 1 -5 5])
       drawnow
       F(k) = getframe;

    end

    movie(F)
    
    end
end