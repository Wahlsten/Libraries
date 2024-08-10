classdef scheme_obj
    % Numerical formulation
    % ut + Au = 0
    properties
        % Numerical scheme
        A;                      % A in Au = F
        F;                      % F in Au = F
        Sigma;                  % Penalty matrices
        Scheme;                 % Scheme matrices:
        
        
        % Numerical solution
        Solution;               % Structure with solutions: u_a, u_n, u_ns, e_L2
        
    end
    methods
        function obj = createAMatrix(obj, grid_obj, sbp_obj, problem_obj)
            % Computes matrix A and F in : AU = F
            
            obj.Scheme.E_0         = sparse(zeros(grid_obj.N(1) + 1));
            obj.Scheme.E_0(1,1)    = 1;
            obj.Scheme.PinvE0Sigma = sparse(sbp_obj.Pinv * obj.Sigma * obj.Scheme.E_0);
            
            force = problem_obj.data.force(grid_obj.t');
            f     = problem_obj.data.f(grid_obj.t');
            
            Dt = sparse(sbp_obj.D1);
            
            obj.A = sparse(Dt - obj.Scheme.PinvE0Sigma);
            obj.F = force - obj.Scheme.PinvE0Sigma * f;
            
        end
        function obj = createPenalties(obj)
            
            obj.Sigma = -1;
            
        end
        function obj = solveUsingBackslash(obj)
            
            obj.Solution.u_n = obj.A\obj.F;
            
        end
        function obj = transformation(obj, problem_obj)
            
            obj.Trans.A = sparse(problem_obj.disc.A);
            
        end
        function plotSolution(obj, grid_obj, problem_obj)
            
            for k = 1:problem_obj.dim
                
                figure
                plot(grid_obj.grid, obj.Solution.u_n(k:problem_obj.dim:end - problem_obj.dim + k, 1))
                
            end
            
        end
        function obj = computeError(obj, grid_obj, sbp_obj, problem_obj)
            
            obj   = computeAnalyticalSolution(obj, grid_obj, problem_obj, grid_obj.maxv(end));
            error = abs(obj.Solution.u_n - obj.Solution.u_a);
            obj.Solution.e_L2 = sqrt(error' * inv(sbp_obj.Pinv) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
                        
            t = grid_obj.t';
            obj.Solution.u_a = problem_obj.data.u_a(t);
                
        end
    end
end