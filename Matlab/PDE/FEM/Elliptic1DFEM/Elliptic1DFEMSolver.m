function [s] = Elliptic1DFEMSolver(problem, acc, data, N, minv, maxv)
%% Comparison

%% Initializing grid
grid          = grid_obj;
grid.minv     = minv;
grid.maxv     = maxv;
grid.Nx        = N;
grid          = createGrid(grid);

%% Initializing problem
p = problem_obj;
p.type.problem = problem;
p.type.data    = data;
p = initialize(p);
p = createMatrices(p, grid);
p = createData(p);
p = createBasis(p);

%% Initializing basis functions
% sbp     = sbp_obj;
% sbp.acc = acc;
% sbp     = createOperators(sbp, grid);

%% Initializing scheme
s = scheme_obj;
s = createScheme(s, grid, p);
s = computeError(s, grid, p);
end