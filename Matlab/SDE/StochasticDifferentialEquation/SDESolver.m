function [s, t] = SDESolver(problem, data, N, minv, maxv, method, R, M)
% Advection Solver

%% Initializing grid§
grid          = grid_obj;
grid.minv     = minv;
grid.maxv     = maxv;
grid.N        = N;
grid          = createGrid(grid);

stc           = stoch_obj;
stc.method    = method; % 'Euler', 'Milstone'
stc.N_mean    = M;
stc.R         = R;
stc           = createGridS(stc, grid);

%% Initializing problem
p = problem_obj;
p.type.problem = problem; % 'Advection1D'
p.type.data    = data; % 'sin'
p = initialize(p);
p = createData(p);

%% Initializing scheme
s = scheme_obj;
tic
s = solveUsingEuler(s, grid, stc, p);
t = toc;
s = computeError(s, grid, p, stc);
% s.plotSolution(grid, stc)
% fprintf(['time elapsed: ', num2str(t), ' sec\n']);
% fprintf(['normed error: ', num2str(s.Solution.error_final), '\n']);
end