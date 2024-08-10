%% Comparison
clear
close all
load 'test25NonSmooth.mat'
NVec    = {0, 0}; %{10:20, 3:2:30};
t       = {zeros(size(NVec{1})), zeros(size(NVec{2}))};
methods = {'NI', 'PC'}; %{'PC', 'NI'};
mean_cmp_normed_t     = {zeros(size(NVec{1})), zeros(size(NVec{2}))};
variance_cmp_normed_t = {zeros(size(NVec{1})), zeros(size(NVec{2}))};

for m = 1:length(methods)
    
    for N = 1:size(NVec{m},2)
        
        disp(NVec{m}(N))
        
        %% Initializing UQ
        uq                         = uq_obj;
        uq.Distribution.minv       = -1;
        uq.Distribution.maxv       =  1;
        uq.GridData.N              = NVec{m}(N); %NVec(m, N);
        uq.GridData.dim            = 1;
        uq.GridGeom.N              = 0;
        uq.GridGeom.dim            = 0;
        uq.method                  = methods{m};
        uq.Distribution.name       = 'Uniform';
        uq.Distribution.parameters = [0, 1];
        uq                         = createGrid(uq);
        
        %% Initializing grid
        grid          = grid_obj;
        grid.minv     = [0, 0];
        grid.maxv     = [1, 1];
        grid.N        = [80 100];
        grid.geometry = 'Unit'; %'Square'; 'Circle Sector'
        grid          = createGrid(grid);
        
        %% Initializing problem
        p = problem_obj;
        p.type.problem = 'Burger'; % 'Burger', 'BurgersSystem', 'BurgersSystemPC7'
        p.type.data    = 'Sin'; % 'BurgersSystemConstant', 'Sin', sincos_smooth', 'sincos_smooth', 'BurgersSystemPC7'
        p.type.uq      = uq.method;
        p = initialize(p);
        p = createMatrices(p, grid, uq);
        p = createData(p);
        
        %% Initializing SBP operators
        sbp     = sbp_obj;
        sbp.acc = 4;
        sbp     = createOperators(sbp, grid);
        sbp     = createDissipationOperator(sbp, grid);
        
        %% Initializing scheme
        s = scheme_obj;
        s = createAMatrix(s, grid, sbp, p);
        s = createScheme(s, grid, sbp, p, uq);
        tic
        s = solveUsingRK4(s, grid, p, uq);
        t{m}(N) = toc;
        s = computeError(s, grid, sbp, p);
        
        %% Results: Statistics
        uq.statistics.mean_a     = mean_a; 
        uq.statistics.variance_a = var_a;
        [mean_cmp_normed_t{m}(N), variance_cmp_normed_t{m}(N), uq] = functionals(uq, grid, s.Solution.u_ns);
%         mean_cmp_normed_t
%         variance_cmp_normed_t
    end
end