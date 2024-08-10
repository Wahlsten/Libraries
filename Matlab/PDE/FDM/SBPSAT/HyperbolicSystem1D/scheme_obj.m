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
        function obj = createAMatrix(obj, grid_obj, sbp_obj, problem_obj)
            % Computes matrix A in : Ut + AU = Pen
            
            E_0               = sparse(zeros(grid_obj.N(1) + 1));
            E_N               = sparse(zeros(grid_obj.N(1) + 1));
            E_0(1, 1)         = 1;
            E_N(end, end)     = 1;
            obj.Matrices.I    = eye(grid_obj.N(1) + 1);
            obj.Matrices.Isys = eye(problem_obj.dim);
            
            obj.Scheme.E_East  = sparse(kron(E_N, obj.Matrices.Isys));
            obj.Scheme.E_West  = sparse(kron(E_0, obj.Matrices.Isys));
            
            obj.Scheme.PinvENSigmaE = sparse(kron(sbp_obj.Pinv, obj.Matrices.Isys) * obj.Sigma.E * obj.Scheme.E_East);
            obj.Scheme.PinvE0SigmaW = sparse(kron(sbp_obj.Pinv, obj.Matrices.Isys) * obj.Sigma.W * obj.Scheme.E_West);
            
            obj.Scheme.PinvE0SigmaWH0  = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_0_p_a));
            obj.Scheme.PinvENSigmaEH0  = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_0_m_a));
            
            obj.Scheme.Amat    = obj.Trans.A;
            
            Am  = obj.Scheme.Amat;
            
            Dx = sparse(kron(sbp_obj.D1, obj.Matrices.Isys));
            
            obj.A = - (Am * Dx) ...
                + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E;
            
            obj.A = sparse(obj.A);
            
        end
        function obj = createScheme(obj, grid_obj, problem_obj)
            
            X = grid_obj.grid;
            
            xlen     = grid_obj.N(1) + 1;
            xzeroVec = grid_obj.minv(1) * ones(1, xlen);
            xoneVec  = grid_obj.maxv(1) * ones(1, xlen);

            Amat     = sparse(kron(obj.Matrices.I, problem_obj.disc.A));

            force = @(t) (   reshape(problem_obj.data.uT_a (X, t), [], 1) ...
                +  Amat    * reshape(problem_obj.data.uX_a (X, t), [], 1));

            obj.cA = @(t, u) (obj.A * u ...
                - obj.Scheme.PinvE0SigmaWH0 * reshape(problem_obj.data.gW (xzeroVec, t), [], 1) ...
                - obj.Scheme.PinvENSigmaEH0 * reshape(problem_obj.data.gE (xoneVec,  t), [], 1) ...
                + force(t));

        end
        function obj = createPenalties(obj)
            
            obj.Sigma.E =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            
        end
        function obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj)

            Am = diag(obj.Trans.A);
            Ix = eye(grid_obj.N(1) + 1);

            obj.Boundary_Operators_a.H_W_0_p_a  = sparse(kron(Ix, blkdiag(eye(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), zeros(problem_obj.Lambda.negDimA))));
            obj.Boundary_Operators_a.H_W_0_m_a  = sparse(kron(Ix, diag(ones(1, problem_obj.Lambda.negDimA), problem_obj.Lambda.posDimA + problem_obj.Lambda.zeroDimA)));

            obj.Boundary_Operators_a.H_E_0_p_a  = sparse(kron(Ix, diag([ones(1, problem_obj.Lambda.posDimA), zeros(1, problem_obj.Lambda.zeroDimA), zeros(1, problem_obj.Lambda.negDimA)])'));
            obj.Boundary_Operators_a.H_E_0_m_a  = sparse(kron(Ix, blkdiag(zeros(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), eye(problem_obj.Lambda.negDimA))));

            obj.Lambda.W_p_a = sparse(kron(Ix, blkdiag(problem_obj.Lambda.Aplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Aplus)))));
            obj.Lambda.E_m_a = sparse(kron(Ix, blkdiag(zeros(problem_obj.dim - length(problem_obj.Lambda.Aminus)), problem_obj.Lambda.Aminus)));

        end
        function obj = createBoundaryOperator(obj, grid_obj, problem_obj)
                            
            obj = transformation(obj, grid_obj, problem_obj);
            obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj);
            Ix   = eye(grid_obj.N(1) + 1);
            Isys = eye(problem_obj.dim);

            H_E_0_p_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_E_0_p_a);
            H_E_0_m_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_E_0_m_a);
            H_W_0_p_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_W_0_p_a);
            H_W_0_m_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_W_0_m_a);

            obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n);
            obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n);
            obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n);
            obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n);

            obj.Lambda.E_m = sparse(obj.Lambda.E_m_a);
            obj.Lambda.W_p = sparse(obj.Lambda.W_p_a);

            obj.Boundary_Operators_n.E = obj.Boundary_Operators_n.H_E_m_n;
            obj.Boundary_Operators_n.W = obj.Boundary_Operators_n.H_W_p_n;
                
        end
        function obj = solveUsingRK4(obj, grid_obj, problem_obj)
            
            X  = grid_obj.grid;
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    
            [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                reshape(problem_obj.data.f(X, 0), [], 1), opts);

            obj.Solution.u_n      = obj.Solution.u_n';

        end
        function obj = transformation(obj, grid_obj, problem_obj)
                        
            Ix = eye(grid_obj.N(1) + 1, grid_obj.N(1) + 1);

            obj.Trans.A    = sparse(kron(Ix, problem_obj.disc.A));
                    
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
            obj.Solution.e_L2 = sqrt(error' * kron(inv(sbp_obj.Pinv), obj.Matrices.Isys) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = grid_obj.grid;

            obj.Solution.u_a = reshape(problem_obj.data.u_a(X, t), [], 1);

        end
    end
end