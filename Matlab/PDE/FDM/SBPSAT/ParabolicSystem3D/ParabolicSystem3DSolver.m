function [s, t] = ParabolicSystem3DSolver(problem, acc, data, N, minv, maxv, geometry)

%% Initializing grid
grid          = grid_obj;
grid.minv     = minv;
grid.maxv     = maxv;
grid.N        = N;
grid.geometry = geometry; %'Square'; 'Circle Sector'
grid          = createGrid(grid);
grid.plotGrid()

%% Initializing problem
p = problem_obj;
p.type.problem = problem; % 'Burgers', 'Advection2D', 'Navier-Stokes', 'Maxwell', 'Navier-Stokes2D'
p.type.data    = data; % 'Burgers', 'sin2pi1D', 'sincos_smooth', 'sin2pi', 'sin', 'sin2DMaxwell', 'NSsinX', 'NSsinY', 'NSConst', 'NSXY'
p = initialize(p);
p = createMatrices(p);
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