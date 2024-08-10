function [s, t] = HyperbolicSystem1DUQSolver(problem, acc, data, N, minv, maxv, geometry)
% Hyperbolic system interface solver
    
    %% Initializing UQ
    tic
    uq                         = uq_obj;
    uq.Distribution.minv       = -1;
    uq.Distribution.maxv       = 1;
    uq.GridData.N              = 200;
    uq.GridData.dim            = 1;
    uq.method                  = {'NI'};
    uq.Distribution.name       = 'Uniform';
    uq.Distribution.parameters = [0, 1];
    uq                         = createGrid(uq);
    
    %% Initializing grid
    grid          = grid_obj;
    grid.minv     = minv;
    grid.maxv     = maxv;
    grid.N        = N;
    grid.geometry = geometry;
    grid          = createGrid(grid);
    
    %% Initializing problem
    p = problem_obj;
    p.type.problem = problem; % 'Advection1D'
    p.type.data    = data; % 'sin'
    p.type.uq      = uq.method;
    p = initialize(p);
    p = createMatrices(p);
    p = createData(p);
    
    %% Initializing SBP operators
    sbp     = sbp_obj;
    sbp.acc = acc;
    sbp     = createOperators(sbp, grid);
    
    %% Initializing scheme
    s = scheme_obj;
    s = createBoundaryOperator(s, grid, p);
    s = createPenalties(s, p, uq);
    s = createAMatrix(s, grid, p, sbp, uq);
    s = createScheme(s, grid, sbp, p, uq);
    %tic
    s = solveUsingRK4(s, grid, p, uq);
    %t = toc;
    %s = computeError(s, grid, sbp, p);
    
    %% Results: Statistics
    uq.statistics.mean_a     = 0; %{0, 0}; %mean_a;
    uq.statistics.variance_a = 0; %{0, 0}; %var_a;
    [mean_cmp_normed_t, variance_cmp_normed_t, uq] ...
        = functionals(uq, grid, s.Solution.u_nsLR);
    toc
end