function [s, t] = AdvectionDiffusionSolverFP(grid_params, problem_params, scheme_params)
% 
% grid_params:
% problem_params:
% scheme_params:

%% Initializing grid
grid_str = CreateGrid(grid_params);

%% Initializing problem
problem_str = initialize(problem_param);

%% Initializing scheme
scheme_str = CreateScheme(grid_str, problem_str, scheme_params);

%% Solve
solution_str = Solve(scheme_str);

end