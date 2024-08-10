function [s, t] = ODESolver(problem, acc, data, N, minv, maxv)
% Advection Solver

%% Initializing grid§
grid          = grid_obj;
grid.minv     = minv;
grid.maxv     = maxv;
grid.N        = N;
grid          = createGrid(grid);

%% Initializing problem
p = problem_obj;
p.type.problem = problem; % 'Advection1D'
p.type.data    = data; % 'sin'
p = initialize(p);
p = createData(p);

%% Initializing SBP operators
sbp     = sbp_obj;
sbp.acc = acc;
sbp     = createOperators(sbp, grid);

%% Initializing scheme
s = scheme_obj;
s = createPenalties(s);
s = createAMatrix(s, grid, sbp, p);
tic
s = solveUsingBackslash(s);
t = toc;
s = computeError(s, grid, sbp, p);
end