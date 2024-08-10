function [s, t] = HyperbolicSystem2DSolver(problem, acc, data, N, minv, maxv, geometry)

%% Initializing grid
grid          = grid_obj;
grid.minv     = minv;
grid.maxv     = maxv;
grid.N        = N;
grid.geometry = geometry;
grid          = createGrid(grid);
grid.plotGrid()

%% Initializing problem
p = problem_obj;
p.type.problem = problem;
p.type.data    = data;
p = initialize(p);
p = createMatrices(p, grid);
p = createData(p);

%% Initializing SBP operators
sbp     = sbp_obj;
sbp.acc = acc;
sbp     = createOperators(sbp, grid);

%% Initializing scheme
s = scheme_obj;
s = createBoundaryOperator(s, grid, sbp, p);
s = createPenalties(s);
s = createAMatrix(s, grid, sbp, p);
s = createScheme(s, grid, p);
tic
s = solveUsingRK4(s, grid, p);
t = toc;
s = computeError(s, grid, sbp, p);