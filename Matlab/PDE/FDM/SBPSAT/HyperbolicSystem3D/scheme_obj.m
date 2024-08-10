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
            
            for k = 1:3 %grid_obj.dim
                
                E_0{k}            = zeros(grid_obj.N(k) + 1);
                E_N{k}            = zeros(grid_obj.N(k) + 1);
                E_0{k}(1, 1)      = 1;
                E_N{k}(end, end)  = 1;
                obj.Matrices.I{k} = eye(grid_obj.N(k) + 1);
                
            end
            
            obj.Matrices.Isys = eye(problem_obj.dim);
            
            Ix   = obj.Matrices.I{1};
            Iy   = obj.Matrices.I{2};
            Iz   = obj.Matrices.I{3};
            Ixyz = kron(kron(Ix, Iy), Iz);
            Isys = obj.Matrices.Isys;
            
            obj.Scheme.E_East  = sparse(kron(kron(kron(E_N{1}, Iy), Iz), Isys));
            obj.Scheme.E_West  = sparse(kron(kron(kron(E_0{1}, Iy), Iz), Isys));
            obj.Scheme.E_North = sparse(kron(kron(kron(Ix, E_N{2}), Iz), Isys));
            obj.Scheme.E_South = sparse(kron(kron(kron(Ix, E_0{2}), Iz), Isys));
            obj.Scheme.E_Front = sparse(kron(kron(kron(Ix, Iy), E_N{2}), Isys));
            obj.Scheme.E_Back  = sparse(kron(kron(kron(Ix, Iy), E_0{2}), Isys));
            
            obj.Scheme.PinvENSigmaE = sparse(kron(kron(kron(sbp_obj.Pinv{1}, Iy), Iz), Isys) * obj.Sigma.E * obj.Scheme.E_East);
            obj.Scheme.PinvE0SigmaW = sparse(kron(kron(kron(sbp_obj.Pinv{1}, Iy), Iz), Isys) * obj.Sigma.W * obj.Scheme.E_West);
            obj.Scheme.PinvENSigmaN = sparse(kron(kron(kron(Ix, sbp_obj.Pinv{2}), Iz), Isys) * obj.Sigma.N * obj.Scheme.E_North);
            obj.Scheme.PinvE0SigmaS = sparse(kron(kron(kron(Ix, sbp_obj.Pinv{2}), Iz), Isys) * obj.Sigma.S * obj.Scheme.E_South);
            obj.Scheme.PinvENSigmaF = sparse(kron(kron(kron(Ix, Iy), sbp_obj.Pinv{3}), Isys) * obj.Sigma.F * obj.Scheme.E_Front);
            obj.Scheme.PinvE0SigmaB = sparse(kron(kron(kron(Ix, Iy), sbp_obj.Pinv{3}), Isys) * obj.Sigma.B * obj.Scheme.E_Back);
            
            obj.Scheme.PinvENSigmaEH0  = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_0_m_a));
            obj.Scheme.PinvE0SigmaWH0  = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_0_p_a));
            obj.Scheme.PinvENSigmaNH0  = sparse(obj.Scheme.PinvENSigmaN * sparse(obj.Boundary_Operators_a.H_N_0_m_a));
            obj.Scheme.PinvE0SigmaSH0  = sparse(obj.Scheme.PinvE0SigmaS * sparse(obj.Boundary_Operators_a.H_S_0_p_a));
            obj.Scheme.PinvENSigmaFH0  = sparse(obj.Scheme.PinvENSigmaF * sparse(obj.Boundary_Operators_a.H_F_0_m_a));
            obj.Scheme.PinvE0SigmaBH0  = sparse(obj.Scheme.PinvE0SigmaB * sparse(obj.Boundary_Operators_a.H_B_0_p_a));

            obj.Scheme.Amat = kron(Ixyz, problem_obj.Lambda.A);
            obj.Scheme.Bmat = kron(Ixyz, problem_obj.Lambda.B);
            obj.Scheme.Cmat = kron(Ixyz, problem_obj.Lambda.C);
            
            Am  = sparse(obj.Scheme.Amat);
            Bm  = sparse(obj.Scheme.Bmat);
            Cm  = sparse(obj.Scheme.Cmat);

            Dx = sparse(kron(kron(kron(sbp_obj.D1{1}, Iy), Iz), Isys));
            Dy = sparse(kron(kron(kron(Ix, sbp_obj.D1{2}), Iz), Isys));
            Dz = sparse(kron(kron(kron(Ix, Iy), sbp_obj.D1{3}), Isys));

            obj.A = - (Am * Dx + Bm * Dy + Cm * Dz) ...
                + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E ...
                + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                + obj.Scheme.PinvENSigmaN * obj.Boundary_Operators_n.N ...
                + obj.Scheme.PinvE0SigmaS * obj.Boundary_Operators_n.S ...
                + obj.Scheme.PinvENSigmaF * obj.Boundary_Operators_n.F ...
                + obj.Scheme.PinvE0SigmaB * obj.Boundary_Operators_n.B;

            obj.A = sparse(obj.A);

        end
        function obj = createScheme(obj, grid_obj, problem_obj)
            
            X      = grid_obj.grid(:, :, :, 1);
            Y      = grid_obj.grid(:, :, :, 2);
            Z      = grid_obj.grid(:, :, :, 3);
            A_a    = obj.Scheme.Amat;
            B_a    = obj.Scheme.Bmat;
            C_a    = obj.Scheme.Cmat;
            dim    = problem_obj.dim;

            force = @(t) (   reshape(reshape(problem_obj.data.uT_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                +  A_a     * reshape(reshape(problem_obj.data.uX_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                +  B_a     * reshape(reshape(problem_obj.data.uY_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                +  C_a     * reshape(reshape(problem_obj.data.uZ_a(X(:), Y(:), Z(:), t), [], dim)', [], 1));

            obj.cA = @(t, u) (obj.A * u ...
                - obj.Scheme.PinvENSigmaEH0  * reshape(reshape(problem_obj.data.gE(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaWH0  * reshape(reshape(problem_obj.data.gW(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvENSigmaNH0  * reshape(reshape(problem_obj.data.gN(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvE0SigmaSH0  * reshape(reshape(problem_obj.data.gS(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaFH0  * reshape(reshape(problem_obj.data.gF(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvE0SigmaBH0  * reshape(reshape(problem_obj.data.gB(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                + force(t));

        end
        function obj = createPenalties(obj)
            
            obj.Sigma.E =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            obj.Sigma.N =  obj.Boundary_Operators_n.H_N_m_n' * obj.Lambda.N_m;
            obj.Sigma.S = -obj.Boundary_Operators_n.H_S_p_n' * obj.Lambda.S_p;
            obj.Sigma.F =  obj.Boundary_Operators_n.H_F_m_n' * obj.Lambda.F_m;
            obj.Sigma.B = -obj.Boundary_Operators_n.H_B_p_n' * obj.Lambda.B_p;

        end
        function obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj)

            Ix   = eye(grid_obj.N(1) + 1);
            Iy   = eye(grid_obj.N(2) + 1);
            Iz   = eye(grid_obj.N(3) + 1);
            Ixyz = kron(kron(Ix, Iy), Iz);

            obj.Boundary_Operators_a.H_E_0_p_a  = sparse(kron(Ixyz, diag([ones(1, problem_obj.Lambda.posDimA), zeros(1, problem_obj.Lambda.zeroDimA), zeros(1, problem_obj.Lambda.negDimA)])'));
            obj.Boundary_Operators_a.H_E_0_m_a  = sparse(kron(Ixyz, blkdiag(zeros(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), eye(problem_obj.Lambda.negDimA))));
            
            obj.Boundary_Operators_a.H_W_0_p_a  = sparse(kron(Ixyz, blkdiag(eye(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), zeros(problem_obj.Lambda.negDimA))));
            obj.Boundary_Operators_a.H_W_0_m_a  = sparse(kron(Ixyz, diag(ones(1, problem_obj.Lambda.negDimA), problem_obj.Lambda.posDimA + problem_obj.Lambda.zeroDimA)));

            obj.Boundary_Operators_a.H_N_0_p_a  = sparse(kron(Ixyz, diag([ones(1, problem_obj.Lambda.posDimB), zeros(1, problem_obj.Lambda.zeroDimB), zeros(1, problem_obj.Lambda.negDimB)])'));
            obj.Boundary_Operators_a.H_N_0_m_a  = sparse(kron(Ixyz, blkdiag(zeros(problem_obj.Lambda.posDimB), zeros(problem_obj.Lambda.zeroDimB), eye(problem_obj.Lambda.negDimB))));

            obj.Boundary_Operators_a.H_S_0_p_a  = sparse(kron(Ixyz, blkdiag(eye(problem_obj.Lambda.posDimB), zeros(problem_obj.Lambda.zeroDimB), zeros(problem_obj.Lambda.negDimB))));
            obj.Boundary_Operators_a.H_S_0_m_a  = sparse(kron(Ixyz, diag(ones(1, problem_obj.Lambda.negDimB), problem_obj.Lambda.posDimB + problem_obj.Lambda.zeroDimB)));
            
            obj.Boundary_Operators_a.H_F_0_p_a  = sparse(kron(Ixyz, diag([ones(1, problem_obj.Lambda.posDimC), zeros(1, problem_obj.Lambda.zeroDimC), zeros(1, problem_obj.Lambda.negDimC)])'));
            obj.Boundary_Operators_a.H_F_0_m_a  = sparse(kron(Ixyz, blkdiag(zeros(problem_obj.Lambda.posDimC), zeros(problem_obj.Lambda.zeroDimC), eye(problem_obj.Lambda.negDimC))));

            obj.Boundary_Operators_a.H_B_0_p_a  = sparse(kron(Ixyz, blkdiag(eye(problem_obj.Lambda.posDimC), zeros(problem_obj.Lambda.zeroDimC), zeros(problem_obj.Lambda.negDimC))));
            obj.Boundary_Operators_a.H_B_0_m_a  = sparse(kron(Ixyz, diag(ones(1, problem_obj.Lambda.negDimC), problem_obj.Lambda.posDimC + problem_obj.Lambda.zeroDimC)));
            
            obj.Lambda.E_m_a = sparse(kron(Ixyz, blkdiag(zeros(problem_obj.dim - length(problem_obj.Lambda.Aminus)), problem_obj.Lambda.Aminus)));
            obj.Lambda.W_p_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.Aplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Aplus)))));
            obj.Lambda.N_m_a = sparse(kron(Ixyz, blkdiag(zeros(problem_obj.dim - length(problem_obj.Lambda.Bminus)), problem_obj.Lambda.Bminus)));
            obj.Lambda.S_p_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.Bplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Bplus)))));
            obj.Lambda.F_m_a = sparse(kron(Ixyz, blkdiag(zeros(problem_obj.dim - length(problem_obj.Lambda.Cminus)), problem_obj.Lambda.Cminus)));
            obj.Lambda.B_p_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.Cplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Cplus)))));
                
        end
        function obj = createBoundaryOperator(obj, grid_obj, sbp_obj, problem_obj)
            
            obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj);

            Isys  = eye(problem_obj.dim);
            Dx    = sbp_obj.D1{1};
            Dy    = sbp_obj.D1{2};
            Dz    = sbp_obj.D1{3};
            Ix    = eye(size(Dx));
            Iy    = eye(size(Dy));
            Iz    = eye(size(Dz));
            Ixyzs = kron(kron(kron(Ix, Iy), Iz), Isys);
            
            H_E_0_p_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_E_0_p_a);
            H_E_0_m_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_E_0_m_a);
            H_W_0_p_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_W_0_p_a);
            H_W_0_m_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_W_0_m_a);
            H_N_0_p_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_N_0_p_a);
            H_N_0_m_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_N_0_m_a);
            H_S_0_p_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_S_0_p_a);
            H_S_0_m_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_S_0_m_a);
            H_F_0_p_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_F_0_p_a);
            H_F_0_m_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_F_0_m_a);
            H_B_0_p_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_B_0_p_a);
            H_B_0_m_n = sparse(Ixyzs * obj.Boundary_Operators_a.H_B_0_m_a);

            obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n);
            obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n);
            obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n);
            obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n);
            obj.Boundary_Operators_n.H_N_p_n = sparse(H_N_0_p_n);
            obj.Boundary_Operators_n.H_N_m_n = sparse(H_N_0_m_n);
            obj.Boundary_Operators_n.H_S_p_n = sparse(H_S_0_p_n);
            obj.Boundary_Operators_n.H_S_m_n = sparse(H_S_0_m_n);
            obj.Boundary_Operators_n.H_F_p_n = sparse(H_F_0_p_n);
            obj.Boundary_Operators_n.H_F_m_n = sparse(H_F_0_m_n);
            obj.Boundary_Operators_n.H_B_p_n = sparse(H_B_0_p_n);
            obj.Boundary_Operators_n.H_B_m_n = sparse(H_B_0_m_n);

            obj.Lambda.E_m = sparse(obj.Lambda.E_m_a);
            obj.Lambda.W_p = sparse(obj.Lambda.W_p_a);
            obj.Lambda.N_m = sparse(obj.Lambda.N_m_a);
            obj.Lambda.S_p = sparse(obj.Lambda.S_p_a);
            obj.Lambda.F_m = sparse(obj.Lambda.F_m_a);
            obj.Lambda.B_p = sparse(obj.Lambda.B_p_a);

            obj.Boundary_Operators_n.E = obj.Boundary_Operators_n.H_E_m_n;
            obj.Boundary_Operators_n.W = obj.Boundary_Operators_n.H_W_p_n;
            obj.Boundary_Operators_n.N = obj.Boundary_Operators_n.H_N_m_n;
            obj.Boundary_Operators_n.S = obj.Boundary_Operators_n.H_S_p_n;
            obj.Boundary_Operators_n.F = obj.Boundary_Operators_n.H_F_m_n;
            obj.Boundary_Operators_n.B = obj.Boundary_Operators_n.H_B_p_n;
            
        end
        function obj = solveUsingRK4(obj, grid_obj, problem_obj)
            
            X  = grid_obj.grid(:, :, :, 1);
            Y  = grid_obj.grid(:, :, :, 2);
            Z  = grid_obj.grid(:, :, :, 3);
            dim = problem_obj.dim;
            
            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
                                            
            [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                reshape(reshape(problem_obj.data.f(X(:), Y(:), Z(:), 0), [], dim)', [], 1), opts);
            
            % reshape(problem_obj.data.f(X(:), Y(:), Z(:), 0)', [], 1)

            obj.Solution.u_n      = obj.Solution.u_n';
                    
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
            obj.Solution.e_L2 = sqrt(error' * kron(kron(kron(inv(sbp_obj.Pinv{1}), inv(sbp_obj.Pinv{2})), inv(sbp_obj.Pinv{2})) , obj.Matrices.Isys) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X   = grid_obj.grid(:,:,:,1);
            Y   = grid_obj.grid(:,:,:,2);
            Z   = grid_obj.grid(:,:,:,3);
            dim = problem_obj.dim;
            
            obj.Solution.u_a = reshape(reshape(problem_obj.data.u_a(X(:), Y(:), Z(:), t), [], dim)', [], 1);
                
        end
    end
end