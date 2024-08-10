function [s, t] = HyperbolicInterfaceSolver(problem, acc, data, N, minv, maxv, geometry)
% Hyperbolic system interface solver

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
s = createPenalties(s);
s = createAMatrix(s, grid, sbp);
s = createScheme(s, grid, p);
tic
s = solveUsingRK4(s, grid, p);
t = toc;
s = computeError(s, grid, sbp, p);

end