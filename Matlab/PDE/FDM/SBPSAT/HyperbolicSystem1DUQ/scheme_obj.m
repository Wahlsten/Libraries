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
        function obj = createAMatrix(obj, grid_obj, problem_obj, sbp_obj, uq_obj)
            % Computes matrix A in : Ut + AU = Pen
            
            N = grid_obj.N(1) + 1;
            
            E_0          = zeros(N);
            E_0(1,1)     = 1;
            E_N          = zeros(N);
            E_N(end,end) = 1;
            
            Isys = speye(problem_obj.dim);
            IsL  = speye(uq_obj.GridData.N(1) + 1);
            
            PinvE0SigmaW  = kron(kron(sbp_obj.Pinv*E_0, Isys), IsL) * kron(obj.Sigma.W, IsL);
            PinvENSigmaE  = kron(kron(sbp_obj.Pinv*E_N, Isys), IsL) * kron(obj.Sigma.E, IsL);
            
            obj.Scheme.Amat = obj.Trans.A;
            
            AL = kron(obj.Scheme.Amat, IsL);
            
            obj.Scheme.Am  = AL;
            
            Dx = sparse(kron(kron(sbp_obj.D1, Isys), IsL));
            
            A11 = - (AL * Dx) ...
                + PinvE0SigmaW * kron(obj.Boundary_Operators_n.W, IsL) ...
                + PinvENSigmaE * kron(obj.Boundary_Operators_n.E, IsL);
                 
            obj.Scheme.B_W = -sparse(PinvE0SigmaW * kron(obj.Boundary_Operators_n.W, IsL));
            obj.Scheme.B_E = -sparse(PinvENSigmaE * kron(obj.Boundary_Operators_n.E, IsL));
            
            obj.A = sparse(A11);
            
            obj.Scheme.Isys = Isys;
            
        end
        function obj = createScheme(obj, grid_obj, sbp_obj, problem_obj, uq_obj)
            
            X = grid_obj.grid;
            S = uq_obj.GridData.grid;
            
            [XX, SS] = meshgrid(X, S);
            
            xzeroVec = grid_obj.minv(1) * ones(1, grid_obj.N(1) + 1);
            xoneVec  = grid_obj.maxv(1) * ones(1, grid_obj.N(1) + 1);
            
            [XX_L, SS_L] = meshgrid(xzeroVec, S);
            [XX_R, SS_R] = meshgrid(xoneVec, S);
            
            obj.cA = @(t,u)(cAfun(obj, grid_obj, uq_obj, problem_obj, XX, SS, XX_L, XX_R, t, u, SS_L, SS_R));
            
        end
        function obj = createPenalties(obj, problem_obj, uq_obj)
            
            IsL = eye(uq_obj.GridData.N + 1);
            
            AL = kron(obj.Trans.A, IsL);
            
            obj.Sigma.E =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            
            obj.Sigma.L =  AL/2; % AL/2
            obj.Sigma.R = -AL/2; % -AR/2
            
        end
        function obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj)
            
            Ix = eye(grid_obj.N(1) + 1);

            obj.Boundary_Operators_a.H_W_0_p_a = sparse(kron(Ix, ...
                blkdiag(eye(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), zeros(problem_obj.Lambda.negDimA))));
            obj.Boundary_Operators_a.H_W_0_m_a = sparse(kron(Ix, ...
                diag(ones(1, problem_obj.Lambda.negDimA), problem_obj.Lambda.posDimA + problem_obj.Lambda.zeroDimA)));

            obj.Boundary_Operators_a.H_E_0_p_a = sparse(kron(Ix, ...
                diag([ones(1, problem_obj.Lambda.posDimA), zeros(1, problem_obj.Lambda.zeroDimA), zeros(1, problem_obj.Lambda.negDimA)])'));
            obj.Boundary_Operators_a.H_E_0_m_a = sparse(kron(Ix, ...
                blkdiag(zeros(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), eye(problem_obj.Lambda.negDimA))));

            obj.Lambda.W_p_a = sparse(kron(Ix, blkdiag(problem_obj.Lambda.Aplus, zeros(problem_obj.dim - length(problem_obj.Lambda.Aplus)))));
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
        function obj = solveUsingRK4(obj, grid_obj, problem_obj, uq_obj)
            
            X = grid_obj.grid;
            S = uq_obj.GridData.grid;
            
            [XX, SS] = meshgrid(X, S);
            
            if strcmp(uq_obj.method{1}, 'PC')
                
                f_L = obj.innerproduct(uq_obj, problem_obj.data.f, XX(1,:), 0);
                f_L = f_L(:);
            
            elseif strcmp(uq_obj.method{1}, 'NI')
                
                f_L = problem_obj.data.f(XX, 0, SS);
                
            end
            opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
            
            [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, f_L, opts);

            obj.Solution.u_nL = zeros(problem_obj.dim * length(X), length(grid_obj.t), uq_obj.GridData.N(1) + 1);

            NL = uq_obj.GridData.N(1)+1;
            
            NXL = (grid_obj.N(1) + 1) * problem_obj.dim;

            for k = 1:NXL

                obj.Solution.u_nsL(:, k, :) = obj.Solution.u_n(:, (k-1) * NL + 1:k*NL);

            end
            
            obj.Solution.u_nsLR = obj.Solution.u_nsL;
            
        end
        function obj = transformation(obj, grid_obj, problem_obj)
            
            Ix = eye(grid_obj.N(1) + 1, grid_obj.N(1) + 1);

            obj.Trans.A = sparse(kron(Ix, problem_obj.disc.A));
            
        end
        function plotSolution(obj, grid_obj, problem_obj)
            
            for k = 1:problem_obj.dim
                
                figure
                plot(grid_obj.grid, obj.Solution.u_n(k:problem_obj.dim:end - problem_obj.dim + k, 1))
                
            end
            
        end
        function obj = computeError(obj, grid_obj, sbp_obj, problem_obj)
            
            obj   = computeAnalyticalSolution(obj, grid_obj, problem_obj, grid_obj.maxv(end));
            error = abs(obj.Solution.u_ns - obj.Solution.u_a);
            obj.Solution.e_L2 = sqrt(error' * inv(blkdiag(kron(sbp_obj.Pinv{1}, obj.Scheme.Isys), kron(sbp_obj.Pinv{2}, obj.Scheme.Isys))) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = [grid_obj.grid{1}, grid_obj.grid{2}];
            
            obj.Solution.u_a = reshape(problem_obj.data.u_a(X, t, 0), [], 1);
                
        end
        function c   = innerproduct(obj, uq_obj, f, x, t)
            
            weight = NodesWeights2(uq_obj, uq_obj.Distribution.name, uq_obj.MaxPol + 1);
            r      = kron(obj.Scheme.Isys, uq_obj.IPfact) * f(x, t, weight(:, 2));
            c      = r(:);
            
        end
        function cA  = cAfun(obj, grid_obj, uq_obj, problem_obj, XX, SS, XX_L, XX_R, t, u, SS_L, SS_R)
            
            if strcmp(uq_obj.method, 'PC')
                
                g_L = obj.innerproduct(uq_obj, problem_obj.data.gW,   XX_L(1, :), t);
                g_R = obj.innerproduct(uq_obj, problem_obj.data.gE,   XX_R(1, :), t);
                ux  = obj.innerproduct(uq_obj, problem_obj.data.uX_a, XX(1, :),   t);
                ut  = obj.innerproduct(uq_obj, problem_obj.data.uT_a, XX(1, :),   t);
                
            elseif strcmp(uq_obj.method, 'NI')
                
                g_L = problem_obj.data.gW(XX_L, t, SS_L);
                g_R = problem_obj.data.gE(XX_R, t, SS_R);
                ux  = problem_obj.data.uX_a(XX, t, SS);
                ut  = problem_obj.data.uT_a(XX, t, SS);
            
            end
            
            g_W = g_L(:);
            g_E = g_R(:);
            
            force = @(t) (ut(:) + obj.Scheme.Am * ux(:));
            
            cAtemp = @(t, u) ( obj.A * u ...
                + obj.Scheme.B_W * g_W ...
                + obj.Scheme.B_E * g_E ...
                + force(t));
            
            cA = cAtemp(t, u);
            
        end
    end
end