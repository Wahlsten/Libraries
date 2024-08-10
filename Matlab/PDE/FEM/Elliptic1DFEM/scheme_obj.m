classdef scheme_obj
    % Numerical formulation
    % ut + Au = 0
    properties
        % Numerical scheme
        A;                      % A in  cA = Au + pen*g
        
        % Numerical solution
        Solution;               % Structure with solutions: u_a, u_n, u_ns, e_L2
        
    end
    methods
        function obj = createScheme(obj, grid_obj, problem_obj)
            % Computes matrix A in : Ut + AU = F
            
            X = grid_obj.grid;
            dx = grid_obj.dx;
            Nx = grid_obj.Nx;
            force = problem_obj.data.force;
            f     = problem_obj.data.f;
            Am    = sparse(zeros(Nx+1));
            M     = sparse(zeros(Nx+1));
            F     = sparse(zeros(Nx+1, 1));
            h1    = problem_obj.data.hat1;
            h2    = problem_obj.data.hat2;
            
            Am(1,1)      = 1;
            F(1)         = f(X(1));
            Am(end, end) = 1;
            F(end)       = f(X(end));
            Am(2, 2)     = 1/dx(1);
            M(2, 2)      = dx(1)/3;
            F(2)         = obj.innerProduct(X(1), X(2), h1, force);
            
            for k = 2:Nx-1
            
                Am(k, k)     = Am(k, k)     + 1/dx(k);
                Am(k, k+1)   = Am(k, k+1)   - 1/dx(k);
                Am(k+1, k)   = Am(k+1, k)   - 1/dx(k);
                Am(k+1, k+1) = Am(k+1, k+1) + 1/dx(k);
                
                M(k, k)     = M(k, k)     + dx(k)/3;
                M(k, k+1)   = M(k, k+1)   + dx(k)/6;
                M(k+1, k)   = M(k+1, k)   + dx(k)/6;
                M(k+1, k+1) = M(k+1, k+1) + dx(k)/3;
                
                F(k)        = F(k)        + obj.innerProduct(X(k), X(k+1), h2, force);
                F(k+1)      = F(k+1)      + obj.innerProduct(X(k), X(k+1), h1, force);

            end
            
            Am(Nx, Nx) = Am(Nx, Nx) + 1/dx(end);
            M(Nx, Nx)  = M(Nx, Nx) + dx(end)/3;
            
            F(Nx) = F(Nx) + obj.innerProduct(X(Nx), X(Nx+1), h2, force);
            obj.A  = problem_obj.disc.A * Am + problem_obj.disc.B * M;
            
            obj.Solution.u_n = obj.A\F;
            plot(X, obj.Solution.u_n, X, problem_obj.data.u_a(X))
    
        end
        function I = innerProduct(~, x1, x2, h, f)
            
            I = integral(@(x)(h(x, x1, x2) .* f(x)), x1, x2);
            
        end
        function plotSolution(obj, grid_obj, problem_obj)
            
            for k = 1:problem_obj.dim
                
                figure
                plot(grid_obj.grid, obj.Solution.u_n(k:problem_obj.dim:end - problem_obj.dim + k, 1))
                
            end
            
        end
        function obj = computeError(obj, grid_obj, problem_obj)
            
            obj   = computeAnalyticalSolution(obj, grid_obj, problem_obj);
            error = abs(obj.Solution.u_n(:,end) - obj.Solution.u_a);
            obj.Solution.e_L2 = norm(error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = grid_obj.grid;
            
            obj.Solution.u_a = reshape(problem_obj.data.u_a(X), [], 1);
                
        end
    end
end