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
            
            E_0            = cell(1, grid_obj.dim + 1);
            E_N            = cell(1, grid_obj.dim + 1);
            obj.Matrices.I = cell(1, grid_obj.dim + 1);
            
            for k = 1:2 %grid_obj.dim
                
                E_0{k}            = zeros(grid_obj.N(k) + 1);
                E_N{k}            = zeros(grid_obj.N(k) + 1);
                E_0{k}(1, 1)      = 1;
                E_N{k}(end, end)  = 1;
                obj.Matrices.I{k} = eye(grid_obj.N(k) + 1);
                
            end
            
            obj.Matrices.Isys = eye(problem_obj.dim);
            
            obj.Scheme.E_East  = sparse(kron(kron(E_N{1}, obj.Matrices.I{2}), obj.Matrices.Isys));
            obj.Scheme.E_West  = sparse(kron(kron(E_0{1}, obj.Matrices.I{2}), obj.Matrices.Isys));
            obj.Scheme.E_North = sparse(kron(kron(obj.Matrices.I{1}, E_N{2}), obj.Matrices.Isys));
            obj.Scheme.E_South = sparse(kron(kron(obj.Matrices.I{1}, E_0{2}), obj.Matrices.Isys));
                            
            obj.Scheme.PinvENSigmaE = sparse(kron(kron(sbp_obj.Pinv{1},   obj.Matrices.I{2}), obj.Matrices.Isys) * obj.Sigma.E * obj.Scheme.E_East);
            obj.Scheme.PinvE0SigmaW = sparse(kron(kron(sbp_obj.Pinv{1},   obj.Matrices.I{2}), obj.Matrices.Isys) * obj.Sigma.W * obj.Scheme.E_West);
            obj.Scheme.PinvENSigmaN = sparse(kron(kron(obj.Matrices.I{1}, sbp_obj.Pinv{2}),   obj.Matrices.Isys) * obj.Sigma.N * obj.Scheme.E_North);
            obj.Scheme.PinvE0SigmaS = sparse(kron(kron(obj.Matrices.I{1}, sbp_obj.Pinv{2}),   obj.Matrices.Isys) * obj.Sigma.S * obj.Scheme.E_South);

            obj.Scheme.PinvE0SigmaWH0  = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_0_p_a));
            obj.Scheme.PinvENSigmaEH0  = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_0_m_a));
            obj.Scheme.PinvE0SigmaSH0  = sparse(obj.Scheme.PinvE0SigmaS * sparse(obj.Boundary_Operators_a.H_S_0_p_a));
            obj.Scheme.PinvENSigmaNH0  = sparse(obj.Scheme.PinvENSigmaN * sparse(obj.Boundary_Operators_a.H_N_0_m_a));

            obj.Scheme.Amat    = obj.Trans.A;
            obj.Scheme.Axmat   = obj.Trans.ADxi;
            obj.Scheme.Bmat    = obj.Trans.B;
            obj.Scheme.Bymat   = obj.Trans.BDyi;

            Am  = obj.Scheme.Amat;
            Amx = obj.Scheme.Axmat;
            Bm  = obj.Scheme.Bmat;
            Bmy = obj.Scheme.Bymat;

            Dx = sparse(kron(kron(sbp_obj.D1{1},     obj.Matrices.I{2}), obj.Matrices.Isys));
            Dy = sparse(kron(kron(obj.Matrices.I{1}, sbp_obj.D1{2}),     obj.Matrices.Isys));

            obj.A = - (Am * Dx + 0*Dx * Am + Amx) ...
                -     (Bm * Dy + 0*Dy * Bm + Bmy) ...
                + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E ...
                + obj.Scheme.PinvE0SigmaS * obj.Boundary_Operators_n.S ...
                + obj.Scheme.PinvENSigmaN * obj.Boundary_Operators_n.N;

            obj.A = sparse(obj.A);

        end
        function obj = createScheme(obj, grid_obj, problem_obj)
            
            X      = grid_obj.gridT(:, :, 1);
            Y      = grid_obj.gridT(:, :, 2);
            A_a    = problem_obj.Trans.A_a;
            B_a    = problem_obj.Trans.B_a;

            force = @(t) (   reshape(problem_obj.data.uT_a (X(:), Y(:), t)', [], 1) ...
                +  A_a     * reshape(problem_obj.data.uX_a (X(:), Y(:), t)', [], 1) ...
                +  B_a     * reshape(problem_obj.data.uY_a (X(:), Y(:), t)', [], 1));

            obj.cA = @(t, u) (obj.A * u ...
                - obj.Scheme.PinvE0SigmaWH0  * reshape(problem_obj.data.gW  (X(:), Y(:), t)', [], 1) ...
                - obj.Scheme.PinvENSigmaEH0  * reshape(problem_obj.data.gE  (X(:), Y(:), t)', [], 1) ...
                - obj.Scheme.PinvE0SigmaSH0  * reshape(problem_obj.data.gS  (X(:), Y(:), t)', [], 1) ...
                - obj.Scheme.PinvENSigmaNH0  * reshape(problem_obj.data.gN  (X(:), Y(:), t)', [], 1) ...
                + force(t));
            
        end
        function obj = createPenalties(obj)
            
            obj.Sigma.E =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            obj.Sigma.N =  obj.Boundary_Operators_n.H_N_m_n' * obj.Lambda.N_m;
            obj.Sigma.S = -obj.Boundary_Operators_n.H_S_p_n' * obj.Lambda.S_p;

        end
        function obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj)

            Ix   = eye(grid_obj.N(1) + 1);
            Iy   = eye(grid_obj.N(2) + 1);

            obj.Boundary_Operators_a.H_W_0_p_a  = sparse(kron(kron(Ix, Iy), blkdiag(eye(problem_obj.Lambda.posDimA),  zeros(problem_obj.Lambda.zeroDimA), zeros(problem_obj.Lambda.negDimA))));
            obj.Boundary_Operators_a.H_W_0_m_a  = sparse(kron(kron(Ix, Iy), diag(ones(1, problem_obj.Lambda.negDimA), problem_obj.Lambda.posDimA + problem_obj.Lambda.zeroDimA)));

            obj.Boundary_Operators_a.H_E_0_p_a  = sparse(kron(kron(Ix, Iy), diag([ones(1, problem_obj.Lambda.posDimA), zeros(1, problem_obj.Lambda.zeroDimA), zeros(1, problem_obj.Lambda.negDimA)])'));
            obj.Boundary_Operators_a.H_E_0_m_a  = sparse(kron(kron(Ix, Iy), blkdiag(zeros(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), eye(problem_obj.Lambda.negDimA))));

            obj.Boundary_Operators_a.H_S_0_p_a  = sparse(kron(kron(Ix, Iy), zeros(problem_obj.dim)));
            obj.Boundary_Operators_a.H_S_0_m_a  = sparse(kron(kron(Ix, Iy), zeros(problem_obj.dim)));
            obj.Boundary_Operators_a.H_N_0_p_a  = sparse(kron(kron(Ix, Iy), zeros(problem_obj.dim)));
            obj.Boundary_Operators_a.H_N_0_m_a  = sparse(kron(kron(Ix, Iy), zeros(problem_obj.dim)));

            obj.Lambda.W_p_a = sparse(kron(kron(Ix, Iy), blkdiag(problem_obj.Lambda.Aplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Aplus)))));
            obj.Lambda.E_m_a = sparse(kron(kron(Ix, Iy), blkdiag(zeros(problem_obj.dim - length(problem_obj.Lambda.Aminus)), problem_obj.Lambda.Aminus)));
            obj.Lambda.S_p_a = sparse(kron(kron(Ix, Iy), blkdiag(problem_obj.Lambda.Bplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Bplus)))));
            obj.Lambda.N_m_a = sparse(kron(kron(Ix, Iy), blkdiag(zeros(problem_obj.dim - length(problem_obj.Lambda.Bminus)), problem_obj.Lambda.Bminus)));
                
        end
        function obj = createBoundaryOperator(obj, grid_obj, sbp_obj, problem_obj)
            
            obj = transformation(obj, grid_obj, sbp_obj, problem_obj);
            obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj);

            Isys = eye(problem_obj.dim);
            Dx   = sbp_obj.D1{1};
            Dy   = sbp_obj.D1{2};
            Ix   = eye(size(Dx));
            Iy   = eye(size(Dy));

            H_E_0_p_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_E_0_p_a);
            H_E_0_m_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_E_0_m_a);
            H_W_0_p_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_W_0_p_a);
            H_W_0_m_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_W_0_m_a);
            H_N_0_p_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_N_0_p_a);
            H_N_0_m_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_N_0_m_a);
            H_S_0_p_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_S_0_p_a);
            H_S_0_m_n = sparse(kron(kron(Ix, Iy), Isys) * obj.Boundary_Operators_a.H_S_0_m_a);

            obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n);
            obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n);
            obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n);
            obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n);
            obj.Boundary_Operators_n.H_N_p_n = sparse(H_N_0_p_n);
            obj.Boundary_Operators_n.H_N_m_n = sparse(H_N_0_m_n);
            obj.Boundary_Operators_n.H_S_p_n = sparse(H_S_0_p_n);
            obj.Boundary_Operators_n.H_S_m_n = sparse(H_S_0_m_n);

            obj.Lambda.E_m = sparse(obj.Lambda.E_m_a);
            obj.Lambda.W_p = sparse(obj.Lambda.W_p_a);
            obj.Lambda.N_m = sparse(obj.Lambda.N_m_a);
            obj.Lambda.S_p = sparse(obj.Lambda.S_p_a);

            obj.Boundary_Operators_n.E = obj.Boundary_Operators_n.H_E_m_n;
            obj.Boundary_Operators_n.W = obj.Boundary_Operators_n.H_W_p_n;
            obj.Boundary_Operators_n.N = obj.Boundary_Operators_n.H_N_m_n;
            obj.Boundary_Operators_n.S = obj.Boundary_Operators_n.H_S_p_n;

        end
        function obj = solveUsingRK4(obj, grid_obj, problem_obj, uq_obj)
            
            X  = grid_obj.gridT(:, :, 1);
            Y  = grid_obj.gridT(:, :, 2);
            
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
                                            
            [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                reshape(problem_obj.data.f(X(:), Y(:), 0)', [], 1), opts);

            obj.Solution.u_n      = obj.Solution.u_n';
                    
        end
        function obj = transformation(obj, grid_obj, sbp_obj, problem_obj, s)
            
            if strcmp(grid_obj.geometry, 'Square')
                
                Ix = eye(grid_obj.N(1) + 1, grid_obj.N(1) + 1);
                Iy = eye(grid_obj.N(2) + 1, grid_obj.N(2) + 1);

                obj.Trans.A    = sparse(kron(kron(Ix, Iy), problem_obj.disc.A));
                obj.Trans.B    = sparse(kron(kron(Ix, Iy), problem_obj.disc.B));
                obj.Trans.ADxi = sparse(kron(kron(Ix, Iy), problem_obj.disc.Ax));
                obj.Trans.BDyi = sparse(kron(kron(Ix, Iy), problem_obj.disc.By));
                
            elseif strcmp(grid_obj.geometry, 'RectangleSine')
                
                D_xi = kron(sbp_obj.D1{1}, eye(size(sbp_obj.D1{2})));
                D_yi = kron(eye(size(sbp_obj.D1{1})), sbp_obj.D1{2});
                A1p  = problem_obj.disc.A;
                A2p  = problem_obj.disc.B;
                NxNy = (grid_obj.N(1) + 1) * (grid_obj.N(2) + 1);
                X = grid_obj.grid(:,:,1);
                Y = grid_obj.grid(:,:,2);
                X = X(:);
                Y = Y(:);
                
                %% New domain
                
                Jdet        =   D_yi * ((D_xi * X) .* Y) - D_xi * ((D_yi * X) .* Y);
                J_xi_x_n    =   D_yi * Y;
                J_xi_y_n    = - D_yi * X;
                J_yi_x_n    = - D_xi * Y;
                J_yi_y_n    =   D_xi * X;
                J_xi_x_xi_n =   D_xi * D_yi * Y;
                J_xi_y_xi_n = - D_xi * D_yi * X;
                J_yi_x_yi_n = - D_yi * D_xi * Y;
                J_yi_y_yi_n =   D_yi * D_xi * X;
                
                obj.Trans.JdetM = spdiags(Jdet,        0, NxNy, NxNy);
                J_xi_x_nM       = spdiags(J_xi_x_n,    0, NxNy, NxNy);
                J_xi_y_nM       = spdiags(J_xi_y_n,    0, NxNy, NxNy);
                J_yi_x_nM       = spdiags(J_yi_x_n,    0, NxNy, NxNy);
                J_yi_y_nM       = spdiags(J_yi_y_n,    0, NxNy, NxNy);
                J_xi_x_xi_nM    = spdiags(J_xi_x_xi_n, 0, NxNy, NxNy);
                J_xi_y_xi_nM    = spdiags(J_xi_y_xi_n, 0, NxNy, NxNy);
                J_yi_x_yi_nM    = spdiags(J_yi_x_yi_n, 0, NxNy, NxNy);
                J_yi_y_yi_nM    = spdiags(J_yi_y_yi_n, 0, NxNy, NxNy);
                
                obj.Trans.A     = sparse(J_xi_x_nM    * A1p + J_xi_y_nM * A2p);
                obj.Trans.B     = sparse(J_yi_x_nM    * A1p + J_yi_y_nM * A2p);
                obj.Trans.ADxi  = sparse(J_xi_x_xi_nM * A1p + J_xi_y_xi_nM * A2p);
                obj.Trans.BDyi  = sparse(J_yi_x_yi_nM * A1p + J_yi_y_yi_nM * A2p);
                
            else
                
                D_xi = kron(sbp_obj.D1{1}, eye(size(sbp_obj.D1{2})));
                D_yi = kron(eye(size(sbp_obj.D1{1})), sbp_obj.D1{2});
                A1p  = problem_obj.disc.A;
                A2p  = problem_obj.disc.B;
                NxNy = (grid_obj.N(1) + 1) * (grid_obj.N(2) + 1);
                ix   = ones(grid_obj.N(1) + 1, 1);
                iy   = ones(grid_obj.N(2) + 1, 1);
                
                Rn  = kron((linspace(grid_obj.Trans.r0(s), grid_obj.Trans.r1(s), grid_obj.N(1) + 1)'), iy );
                FIn = kron(ix, (linspace(grid_obj.Trans.fi0(s), grid_obj.Trans.fi1(s), grid_obj.N(2) + 1)'));
                
                R   = Rn  - 0.05 * s * sin(2 * pi * (FIn - grid_obj.Trans.fi0(s))./(grid_obj.Trans.fi1(s) - grid_obj.Trans.fi0(s)));
                FI  = FIn + 0.05 * s * sin(2 * pi * (Rn  - grid_obj.Trans.r0 (s))./(grid_obj.Trans.r1 (s) - grid_obj.Trans.r0 (s)));
                
                %% New domain
                
                X = R .* cos(FI);
                Y = R .* sin(FI);
                
                Jdet        =   D_yi * ((D_xi * X) .* Y) - D_xi * ((D_yi * X) .* Y);
                J_xi_x_n    =   D_yi * Y;
                J_xi_y_n    = - D_yi * X;
                J_yi_x_n    = - D_xi * Y;
                J_yi_y_n    =   D_xi * X;
                J_xi_x_xi_n =   D_xi * D_yi * Y;
                J_xi_y_xi_n = - D_xi * D_yi * X;
                J_yi_x_yi_n = - D_yi * D_xi * Y;
                J_yi_y_yi_n =   D_yi * D_xi * X;
                
                obj.Trans.JdetM = spdiags(Jdet,        0, NxNy, NxNy);
                J_xi_x_nM       = spdiags(J_xi_x_n,    0, NxNy, NxNy);
                J_xi_y_nM       = spdiags(J_xi_y_n,    0, NxNy, NxNy);
                J_yi_x_nM       = spdiags(J_yi_x_n,    0, NxNy, NxNy);
                J_yi_y_nM       = spdiags(J_yi_y_n,    0, NxNy, NxNy);
                J_xi_x_xi_nM    = spdiags(J_xi_x_xi_n, 0, NxNy, NxNy);
                J_xi_y_xi_nM    = spdiags(J_xi_y_xi_n, 0, NxNy, NxNy);
                J_yi_x_yi_nM    = spdiags(J_yi_x_yi_n, 0, NxNy, NxNy);
                J_yi_y_yi_nM    = spdiags(J_yi_y_yi_n, 0, NxNy, NxNy);
                
                obj.Trans.A     = sparse(J_xi_x_nM    * A1p + J_xi_y_nM * A2p);
                obj.Trans.B     = sparse(J_yi_x_nM    * A1p + J_yi_y_nM * A2p);
                obj.Trans.ADxi  = sparse(J_xi_x_xi_nM * A1p + J_xi_y_xi_nM * A2p);
                obj.Trans.BDyi  = sparse(J_yi_x_yi_nM * A1p + J_yi_y_yi_nM * A2p);
                
            end
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
            obj.Solution.e_L2 = sqrt(error' * kron(kron(inv(sbp_obj.Pinv{1}), inv(sbp_obj.Pinv{2})), obj.Matrices.Isys) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = grid_obj.gridT(:,:,1);
            Y = grid_obj.gridT(:,:,2);
            
            if grid_obj.dim == 1
                
                obj.Solution.u_a = reshape(problem_obj.data.u_a(X, t, 0), [], 1);
                
            elseif grid_obj.dim == 2
                
                obj.Solution.u_a = reshape(problem_obj.data.u_a(X(:), Y(:), t, 0)', [], 1);
                
            end
        end
    end
end