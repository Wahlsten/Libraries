function [s, t] = NonLinearSystemParabolic1DSolver(problem, acc, data, N, minv, maxv, geometry)

variance_v = zeros(40, 1);
mean_v     = zeros(40, 1);

for k = 1:40
    
    disp(k)
    
    %% Initializing UQ
    uq                         = uq_obj;
    uq.Distribution.minv       = -1;
    uq.Distribution.maxv       =  1;
    
    if k == 0
        uq.GridData.N          = 25;
    else
        uq.GridData.N          = 2; %2*k; %k - 1;
    end
    
    uq.GridData.dim            = 1;
    uq.GridGeom.N              = 0;
    uq.GridGeom.dim            = 0;
    uq.method                  = {'PC'};
    uq.Distribution.name       = 'Uniform';
    uq.Distribution.parameters = [0, 1];
    uq                         = createGrid(uq);
    
    %% Initializing grid
    grid          = grid_obj;
    grid.minv     = minv;
    grid.maxv     = maxv;
    grid.N        = N;
    grid.geometry = geometry; %'Unit', 'Square'; 'Circle Sector'
    grid          = createGrid(grid);
    
    %% Initializing problem
    p              = problem_obj;
    p.type.problem = problem; % 'Burger001', 'Burger', 'BurgersSystem'
    p.type.data    = data; % 'BurgersSystemConstant', 'Sin', sincos_smooth', 'sincos_nonsmooth'
    p.type.uq      = uq.method;
    p              = initialize(p);
    p              = createMatrices(p, grid, uq);
    p              = createData(p);
    
    %% Initializing SBP operators
    sbp     = sbp_obj;
    sbp.acc = acc;
    sbp     = createOperators(sbp, grid);
    sbp     = createDissipationOperator(sbp, grid);
    
    %% Initializing scheme
    s = scheme_obj;
    s = createAMatrix(s, grid, sbp, p);
    s = createScheme(s, grid, sbp, p, uq);
    tic
    s = solveUsingRK4(s, grid, p, uq);
    
    if k == 0
        toc
    else
        t(k) = toc;
    end
    %s = computeError(s, grid, sbp, p);
    
    %% Results: Statistics
    %load('PC_25_slow', 'mean_a', 'var_a')
    %load('PC_25_faster', 'mean_a', 'var_a')
    if k == 0
        
        uq.statistics.mean_a     = 0; %mean_a;
        uq.statistics.variance_a = 0; %var_a;
        [~, ~, uq] = functionals(uq, grid, s.Solution.u_ns);
        
        mean_a = uq.statistics.mean;
        var_a = uq.statistics.variance;
        
        save('PC_25_slow', 'mean_a', 'var_a')
        
    else
        
    load('PC_25_slow', 'mean_a', 'var_a')
    uq.statistics.mean_a     = mean_a; %mean_a;
    uq.statistics.variance_a = var_a; %var_a;
    [mean_v(k), variance_v(k), uq] ...
        = functionals(uq, grid, s.Solution.u_ns);
    format long
    variance_v
    
    end
end

if strcmp(uq.method, 'PC')
    
    N_PC = 1:21;
    t_PC = t;
    variance_PC = variance_v;
    save('PC_slow_0_20', 'N_PC', 't_PC', 'variance_PC')
    
elseif strcmp(uq.method, 'NI')
   
    N_NI = 2:2:80;
    t_NI = t;
    variance_NI = variance_v;
    save('NI_slow_2_80', 'N_NI', 't_NI', 'variance_NI')
    
end
end