function [s, t] = AdvectionScalar3DSolver(problem, acc, data, N, minv, maxv, geometry)

tic
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
p.type.problem = problem; % 'Scalar'
p.type.data    = data; % 'Sin', 'constant'
p = initialize(p);
p = createMatrices(p);
p = createData(p);

%% Initializing SBP operators
sbp     = sbp_obj;
sbp.acc = acc;
sbp     = createOperators(sbp, grid);

%% Initializing scheme
s = scheme_obj;
s = createBoundaryOperator(s, p, grid);
s = createPenalties(s);
s = createAMatrix(s, grid, sbp, p);
s = createScheme(s, grid, p);
s = solveUsingRK4(s, grid, p);
t = toc;
s = computeError(s, grid, sbp, p);
%end