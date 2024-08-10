classdef scheme_obj
    % Numerical formulation
    % ut + Au = 0
    properties
        % Numerical scheme
        A;                      % A in  cA = Au + pen*g
        cA;                     % cA in cA = Au + pen*g
        Sigma;                  % Penalty matrices
        Matrices;               % Matrices: I, Isys
        Lambda;                 % Eigenvalues: W_p, E_m
        Boundary_Operators_n;   % Discrete boundary operators: L_W, L_E
        Boundary_Operators_a;   % Continuous boundary operators:
        Scheme;                 % Scheme matrices:
        Trans;                  % Transformation: A, B, C, D, JdetM
        
        % Numerical solution
        Solution;               % Structure with solutions: u_a, u_n, u_ns, e_L2
        
    end
    methods
        function obj = createAMatrix(obj, grid_obj, sbp_obj)
            % Computes matrix A in : Ut + AU = Pen
            
            N_L = grid_obj.N(1) + 1;
            N_R = grid_obj.N(2) + 1;
            
            obj.Scheme.E_0{1}          = zeros(N_L);
            obj.Scheme.E_0{1}(1,1)     = 1;
            obj.Scheme.E_N{1}          = zeros(N_L);
            obj.Scheme.E_N{1}(end,end) = 1;
            
            obj.Scheme.E_0{2}          = zeros(N_R);
            obj.Scheme.E_0{2}(1,1)     = 1;
            obj.Scheme.E_N{2}          = zeros(N_R);
            obj.Scheme.E_N{2}(end,end) = 1;
            
            obj.Scheme.E_I{1}                   = zeros(N_L, N_L + N_R);
            obj.Scheme.E_I{1}(N_L, N_L:N_L + 1) = [1, -1];
            obj.Scheme.E_I{2}                   = zeros(N_R, N_L + N_R);
            obj.Scheme.E_I{2}(1,   N_L:N_L + 1) = [-1, 1];
            
            obj.Scheme.PinvE0SigmaW = sparse(blkdiag(sbp_obj.Pinv{1} * obj.Sigma.W * obj.Scheme.E_0{1}, sparse(N_R, N_R)));
            obj.Scheme.PinvENSigmaE = sparse(blkdiag(sparse(N_L, N_L), sbp_obj.Pinv{2} * obj.Sigma.E * obj.Scheme.E_N{2}));
            
            PinvEISigmaWI = sparse(sbp_obj.Pinv{1} * obj.Sigma.WI * obj.Scheme.E_I{1});
            PinvEISigmaEI = sparse(sbp_obj.Pinv{2} * obj.Sigma.EI * obj.Scheme.E_I{2});
            
            obj.Scheme.Interface = [PinvEISigmaWI; PinvEISigmaEI];
            
            obj.Scheme.PinvE0SigmaWH0 = sparse(obj.Scheme.PinvE0SigmaW ...
                * blkdiag(sparse(obj.Boundary_Operators_a.H_W_0_p_a), sparse(N_R, N_R)));
            
            obj.Scheme.PinvENSigmaEH0 = sparse(obj.Scheme.PinvENSigmaE ...
                * blkdiag(sparse(N_L, N_L), sparse(obj.Boundary_Operators_a.H_E_0_m_a)));
            
            obj.Scheme.Amat = obj.Trans.A;
            obj.Scheme.Bmat = obj.Trans.B;
            
            Am  = blkdiag(obj.Scheme.Amat, obj.Scheme.Bmat);
            
            Dx{1} = sparse(sbp_obj.D1{1});
            Dx{2} = sparse(sbp_obj.D1{2});
            
            D = blkdiag(Dx{1}, Dx{2});
            
            obj.A = - (Am * D) ...
                + obj.Scheme.PinvE0SigmaW * blkdiag(obj.Boundary_Operators_n.W, sparse(N_R, N_R)) ...
                + obj.Scheme.PinvENSigmaE * blkdiag(sparse(N_L, N_L), obj.Boundary_Operators_n.E) ...
                + obj.Scheme.Interface;
            
            obj.A = sparse(obj.A);
            
        end
        function obj = createScheme(obj, grid_obj, problem_obj)
            
            X = [grid_obj.grid{1}, grid_obj.grid{2}];
            
            xlen     = grid_obj.N(1) + 1 + grid_obj.N(2) + 1;
            xzeroVec = grid_obj.minv(1) * ones(1, xlen);
            xoneVec  = grid_obj.maxv(2) * ones(1, xlen);
            Amat     = sparse(problem_obj.disc.A);
            Bmat     = sparse(problem_obj.disc.B);
            AM       = blkdiag(Amat, Bmat);
            
            force = @(t) (reshape(problem_obj.data.uT_a (X, t), [], 1) ...
                +  AM * reshape(problem_obj.data.uX_a (X, t), [], 1));            
            obj.cA = @(t, u) (obj.A * u ...
                - obj.Scheme.PinvE0SigmaWH0  * reshape(problem_obj.data.gW (xzeroVec, t), [], 1) ...
                - obj.Scheme.PinvENSigmaEH0  * reshape(problem_obj.data.gE (xoneVec,  t), [], 1) ...
                + force(t));
            
        end
        function obj = createPenalties(obj)
            
            obj.Sigma.E  =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W  = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            obj.Sigma.EI = -obj.Trans.B(end)/2;
            obj.Sigma.WI =  obj.Trans.A(end)/2;

        end
        function obj = createContinuousBoundaryOperators(obj)
            
            Am = diag(obj.Trans.A);
            
            obj.Boundary_Operators_a.H_E_0_p_a  = sparse(diag( (abs(Am) < 1e-12) ...
                .* abs(Am) + (Am >  1e-12) .* Am));
            obj.Boundary_Operators_a.H_E_0_m_a  = sparse(diag( (abs(Am) < 1e-12) ...
                .* abs(Am) + (Am < -1e-12) .* Am));
            obj.Boundary_Operators_a.H_W_0_p_a  = sparse(diag( (abs(Am) < 1e-12) ...
                .* abs(Am) + (Am >  1e-12) .* Am));
            obj.Boundary_Operators_a.H_W_0_m_a  = sparse(diag( (abs(Am) < 1e-12) ...
                .* abs(Am) + (Am < -1e-12) .* Am));
            
            obj.Lambda.E_m_a = sparse(diag(-pinv(diag(abs(full(Am)))) * (abs(Am) > 1e-12)));
            obj.Lambda.W_p_a = sparse(diag( pinv(diag(abs(full(Am)))) * (abs(Am) > 1e-12)));
            
        end
        function obj = createBoundaryOperator(obj, problem_obj)
            
            obj = transformation(obj, problem_obj);
            obj = createContinuousBoundaryOperators(obj);
            
            H_E_0_p_n = sparse(obj.Boundary_Operators_a.H_E_0_p_a);
            H_E_0_m_n = sparse(obj.Boundary_Operators_a.H_E_0_m_a);
            H_W_0_p_n = sparse(obj.Boundary_Operators_a.H_W_0_p_a);
            H_W_0_m_n = sparse(obj.Boundary_Operators_a.H_W_0_m_a);
            
            obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n);
            obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n);
            obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n);
            obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n);
            
            obj.Lambda.E_m = sparse(diag(-pinv(abs(full(diag(obj.Trans.A)))) ...
                * (abs(diag(obj.Trans.A)) > 1e-12)));
            obj.Lambda.W_p = sparse(diag( pinv(abs(full(diag(obj.Trans.A)))) ...
                * (abs(diag(obj.Trans.A)) > 1e-12)));
            
            obj.Boundary_Operators_n.E = obj.Boundary_Operators_n.H_E_m_n;
            obj.Boundary_Operators_n.W = obj.Boundary_Operators_n.H_W_p_n;
            
        end
        function obj = solveUsingRK4(obj, grid_obj, problem_obj)
            
            X  = [grid_obj.grid{1}, grid_obj.grid{2}];
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
            
            
            [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                reshape(problem_obj.data.f(X, 0), [], 1), opts);
            
            obj.Solution.u_n      = obj.Solution.u_n';
            
        end
        function obj = transformation(obj, problem_obj)
            
            obj.Trans.A = sparse(problem_obj.disc.A);
            obj.Trans.B = sparse(problem_obj.disc.B);
            
        end
        function plotSolution(obj, grid_obj, problem_obj)
            
            for k = 1:problem_obj.dim
                
                figure
                plot(grid_obj.grid, obj.Solution.u_n(k:problem_obj.dim:end - problem_obj.dim + k, 1))
                
            end
            
        end
        function obj = computeError(obj, grid_obj, sbp_obj, problem_obj)
            
            obj   = computeAnalyticalSolution(obj, grid_obj, problem_obj, grid_obj.maxv(end));
            error = abs(obj.Solution.u_n(:,end) - obj.Solution.u_a);
            obj.Solution.e_L2 = sqrt(error' * inv(blkdiag(sbp_obj.Pinv{1}, sbp_obj.Pinv{2})) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = [grid_obj.grid{1}, grid_obj.grid{2}];
            
            obj.Solution.u_a = reshape(problem_obj.data.u_a(X, t), [], 1);
                
        end
    end
end