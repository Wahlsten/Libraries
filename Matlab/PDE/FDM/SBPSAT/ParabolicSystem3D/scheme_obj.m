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
            
            obj.Scheme.PinvE0SigmaWHDx = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_Dx_p_a));
            obj.Scheme.PinvE0SigmaWHDy = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_Dy_p_a));
            obj.Scheme.PinvE0SigmaWHDz = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_Dz_p_a));
            obj.Scheme.PinvENSigmaEHDx = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_Dx_m_a));
            obj.Scheme.PinvENSigmaEHDy = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_Dy_m_a));
            obj.Scheme.PinvENSigmaEHDz = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_Dz_m_a));
            obj.Scheme.PinvE0SigmaSHDx = sparse(obj.Scheme.PinvE0SigmaS * sparse(obj.Boundary_Operators_a.H_S_Dx_p_a));
            obj.Scheme.PinvE0SigmaSHDy = sparse(obj.Scheme.PinvE0SigmaS * sparse(obj.Boundary_Operators_a.H_S_Dy_p_a));
            obj.Scheme.PinvE0SigmaSHDz = sparse(obj.Scheme.PinvE0SigmaS * sparse(obj.Boundary_Operators_a.H_S_Dz_p_a));
            obj.Scheme.PinvENSigmaNHDx = sparse(obj.Scheme.PinvENSigmaN * sparse(obj.Boundary_Operators_a.H_N_Dx_m_a));
            obj.Scheme.PinvENSigmaNHDy = sparse(obj.Scheme.PinvENSigmaN * sparse(obj.Boundary_Operators_a.H_N_Dy_m_a));
            obj.Scheme.PinvENSigmaNHDz = sparse(obj.Scheme.PinvENSigmaN * sparse(obj.Boundary_Operators_a.H_N_Dz_m_a));
            obj.Scheme.PinvE0SigmaBHDx = sparse(obj.Scheme.PinvE0SigmaB * sparse(obj.Boundary_Operators_a.H_B_Dx_p_a));
            obj.Scheme.PinvE0SigmaBHDy = sparse(obj.Scheme.PinvE0SigmaB * sparse(obj.Boundary_Operators_a.H_B_Dy_p_a));
            obj.Scheme.PinvE0SigmaBHDz = sparse(obj.Scheme.PinvE0SigmaB * sparse(obj.Boundary_Operators_a.H_B_Dz_p_a));
            obj.Scheme.PinvENSigmaFHDx = sparse(obj.Scheme.PinvENSigmaF * sparse(obj.Boundary_Operators_a.H_F_Dx_m_a));
            obj.Scheme.PinvENSigmaFHDy = sparse(obj.Scheme.PinvENSigmaF * sparse(obj.Boundary_Operators_a.H_F_Dy_m_a));
            obj.Scheme.PinvENSigmaFHDz = sparse(obj.Scheme.PinvENSigmaF * sparse(obj.Boundary_Operators_a.H_F_Dz_m_a));
            
            obj.Scheme.Amat = kron(Ixyz, problem_obj.disc.A);
            obj.Scheme.Bmat = kron(Ixyz, problem_obj.disc.B);
            obj.Scheme.Cmat = kron(Ixyz, problem_obj.disc.C);
            
            obj.Scheme.D11mat = kron(Ixyz, problem_obj.disc.D11);
            obj.Scheme.D22mat = kron(Ixyz, problem_obj.disc.D22);
            obj.Scheme.D33mat = kron(Ixyz, problem_obj.disc.D33);
            
            obj.Scheme.D12mat = kron(Ixyz, problem_obj.disc.D12);
            obj.Scheme.D13mat = kron(Ixyz, problem_obj.disc.D13);
            obj.Scheme.D21mat = kron(Ixyz, problem_obj.disc.D21);
            obj.Scheme.D23mat = kron(Ixyz, problem_obj.disc.D23);
            obj.Scheme.D31mat = kron(Ixyz, problem_obj.disc.D31);
            obj.Scheme.D32mat = kron(Ixyz, problem_obj.disc.D32);
            
            obj.Scheme.epsilon = problem_obj.epsilon;

            Am   = obj.Scheme.Amat;
            Bm   = obj.Scheme.Bmat;
            Cm   = obj.Scheme.Cmat;
            D11m = obj.Scheme.D11mat;
            D22m = obj.Scheme.D22mat;
            D33m = obj.Scheme.D33mat;
            D12m = obj.Scheme.D12mat;
            D13m = obj.Scheme.D13mat;
            D21m = obj.Scheme.D21mat;
            D23m = obj.Scheme.D23mat;
            D31m = obj.Scheme.D31mat;
            D32m = obj.Scheme.D32mat;
            em   = obj.Scheme.epsilon;

            obj.Scheme.epsD11mat = sparse(problem_obj.epsilon * obj.Scheme.D11mat);
            obj.Scheme.epsD22mat = sparse(problem_obj.epsilon * obj.Scheme.D22mat);
            obj.Scheme.epsD33mat = sparse(problem_obj.epsilon * obj.Scheme.D33mat);
            
            obj.Scheme.epsD12mat = sparse(problem_obj.epsilon * obj.Scheme.D12mat);
            obj.Scheme.epsD13mat = sparse(problem_obj.epsilon * obj.Scheme.D13mat);
            obj.Scheme.epsD21mat = sparse(problem_obj.epsilon * obj.Scheme.D21mat);
            obj.Scheme.epsD23mat = sparse(problem_obj.epsilon * obj.Scheme.D23mat);
            obj.Scheme.epsD31mat = sparse(problem_obj.epsilon * obj.Scheme.D31mat);
            obj.Scheme.epsD32mat = sparse(problem_obj.epsilon * obj.Scheme.D32mat);

            Dx = sparse(kron(kron(kron(sbp_obj.D1{1}, obj.Matrices.I{2}), obj.Matrices.I{3}), obj.Matrices.Isys));
            Dy = sparse(kron(kron(kron(obj.Matrices.I{1}, sbp_obj.D1{2}), obj.Matrices.I{3}), obj.Matrices.Isys));
            Dz = sparse(kron(kron(kron(obj.Matrices.I{1}, obj.Matrices.I{2}), sbp_obj.D1{3}), obj.Matrices.Isys));

            obj.A = - (Am * Dx + Bm * Dy + Cm * Dz) ...
                +      D11m * em * Dx * Dx          ...
                +      D22m * em * Dy * Dy          ...
                +      D33m * em * Dz * Dz          ...
                +      D12m * em * Dx * Dy          ...
                +      D13m * em * Dx * Dz          ...
                +      D21m * em * Dy * Dx          ...
                +      D23m * em * Dy * Dz          ...
                +      D31m * em * Dz * Dx          ...
                +      D32m * em * Dz * Dy          ...
                + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E ...
                + obj.Scheme.PinvE0SigmaS * obj.Boundary_Operators_n.S ...
                + obj.Scheme.PinvENSigmaN * obj.Boundary_Operators_n.N ...
                + obj.Scheme.PinvE0SigmaB * obj.Boundary_Operators_n.B ...
                + obj.Scheme.PinvENSigmaF * obj.Boundary_Operators_n.F;

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
            epsD11_a = problem_obj.epsilon * obj.Scheme.D11mat;
            epsD22_a = problem_obj.epsilon * obj.Scheme.D22mat;
            epsD33_a = problem_obj.epsilon * obj.Scheme.D33mat;
            epsD12_a = problem_obj.epsilon * obj.Scheme.D12mat;
            epsD13_a = problem_obj.epsilon * obj.Scheme.D13mat;
            epsD21_a = problem_obj.epsilon * obj.Scheme.D21mat;
            epsD23_a = problem_obj.epsilon * obj.Scheme.D23mat;
            epsD31_a = problem_obj.epsilon * obj.Scheme.D31mat;
            epsD32_a = problem_obj.epsilon * obj.Scheme.D32mat;

            force = @(t) (    reshape(reshape(problem_obj.data.uT_a (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                +  A_a      * reshape(reshape(problem_obj.data.uX_a (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                +  B_a      * reshape(reshape(problem_obj.data.uY_a (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                +  C_a      * reshape(reshape(problem_obj.data.uZ_a (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD11_a * reshape(reshape(problem_obj.data.uXX_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD22_a * reshape(reshape(problem_obj.data.uYY_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD33_a * reshape(reshape(problem_obj.data.uZZ_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD12_a * reshape(reshape(problem_obj.data.uXY_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD13_a * reshape(reshape(problem_obj.data.uXZ_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD21_a * reshape(reshape(problem_obj.data.uYX_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD23_a * reshape(reshape(problem_obj.data.uYZ_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD31_a * reshape(reshape(problem_obj.data.uZX_a(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                -  epsD32_a * reshape(reshape(problem_obj.data.uZY_a(X(:), Y(:), Z(:), t), [], dim)', [], 1));

            obj.cA = @(t, u) (obj.A * u ...        
                - obj.Scheme.PinvENSigmaEH0  * reshape(reshape(problem_obj.data.gE  (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvENSigmaEHDx * reshape(reshape(problem_obj.data.gEDx(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvENSigmaEHDy * reshape(reshape(problem_obj.data.gEDy(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvENSigmaEHDz * reshape(reshape(problem_obj.data.gEDz(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaWH0  * reshape(reshape(problem_obj.data.gW  (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaWHDx * reshape(reshape(problem_obj.data.gWDx(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaWHDy * reshape(reshape(problem_obj.data.gWDy(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaWHDz * reshape(reshape(problem_obj.data.gWDz(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvENSigmaNH0  * reshape(reshape(problem_obj.data.gN  (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaNHDx * reshape(reshape(problem_obj.data.gNDx(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaNHDy * reshape(reshape(problem_obj.data.gNDy(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaNHDz * reshape(reshape(problem_obj.data.gNDz(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvE0SigmaSH0  * reshape(reshape(problem_obj.data.gS  (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvE0SigmaSHDx * reshape(reshape(problem_obj.data.gSDx(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvE0SigmaSHDy * reshape(reshape(problem_obj.data.gSDy(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvE0SigmaSHDz * reshape(reshape(problem_obj.data.gSDz(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvENSigmaFH0  * reshape(reshape(problem_obj.data.gF  (X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaFHDx * reshape(reshape(problem_obj.data.gFDx(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaFHDy * reshape(reshape(problem_obj.data.gFDy(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...        
                - obj.Scheme.PinvENSigmaFHDz * reshape(reshape(problem_obj.data.gFDz(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaBH0  * reshape(reshape(problem_obj.data.gB  (X(:), Y(:), Z(:), t), [], dim)', [], 1) ... 
                - obj.Scheme.PinvE0SigmaBHDx * reshape(reshape(problem_obj.data.gBDx(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaBHDy * reshape(reshape(problem_obj.data.gBDy(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...
                - obj.Scheme.PinvE0SigmaBHDz * reshape(reshape(problem_obj.data.gBDz(X(:), Y(:), Z(:), t), [], dim)', [], 1) ...       
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
            em   = problem_obj.epsilon;
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
            
            
            
            obj.Boundary_Operators_a.H_W_Dx_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Aplus)  * problem_obj.disc.D11(1:problem_obj.Lambda.posDimA            , :); problem_obj.X.D11plus']));
            obj.Boundary_Operators_a.H_W_Dx_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Aminus) * problem_obj.disc.D11(end - problem_obj.Lambda.negDimA + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D11minus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim); problem_obj.X.D11minus']));

            obj.Boundary_Operators_a.H_E_Dx_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D11plus' ;-inv(problem_obj.Lambda.Aplus) * problem_obj.disc.D11(1:problem_obj.Lambda.posDimA            , :)]));
            obj.Boundary_Operators_a.H_E_Dx_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D11minus'; zeros(problem_obj.dim - size(problem_obj.X.D11minus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Aminus) * problem_obj.disc.D11(end - problem_obj.Lambda.negDimA + 1:end, :)]));

            obj.Boundary_Operators_a.H_S_Dx_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Bplus)  * problem_obj.disc.D12(1:problem_obj.Lambda.posDimB            , :); problem_obj.X.D12plus']));
            obj.Boundary_Operators_a.H_S_Dx_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Bminus) * problem_obj.disc.D12(end - problem_obj.Lambda.negDimB + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D12minus, 2) - length(problem_obj.Lambda.Bminus), problem_obj.dim); problem_obj.X.D12minus']));

            obj.Boundary_Operators_a.H_N_Dx_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D12plus' ;-inv(problem_obj.Lambda.Bplus) * problem_obj.disc.D12(1:problem_obj.Lambda.posDimB            , :)]));
            obj.Boundary_Operators_a.H_N_Dx_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D12minus'; zeros(problem_obj.dim - size(problem_obj.X.D12minus, 2) - length(problem_obj.Lambda.Bminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Bminus) * problem_obj.disc.D12(end - problem_obj.Lambda.negDimB + 1:end, :)]));

            obj.Boundary_Operators_a.H_B_Dx_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Cplus)  * problem_obj.disc.D13(1:problem_obj.Lambda.posDimC            , :); problem_obj.X.D13plus']));
            obj.Boundary_Operators_a.H_B_Dx_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Cminus) * problem_obj.disc.D13(end - problem_obj.Lambda.negDimC + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D13minus, 2) - length(problem_obj.Lambda.Cminus), problem_obj.dim); problem_obj.X.D13minus']));

            obj.Boundary_Operators_a.H_F_Dx_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D13plus' ;-inv(problem_obj.Lambda.Cplus) * problem_obj.disc.D13(1:problem_obj.Lambda.posDimC            , :)]));
            obj.Boundary_Operators_a.H_F_Dx_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D13minus'; zeros(problem_obj.dim - size(problem_obj.X.D13minus, 2) - length(problem_obj.Lambda.Cminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Cminus) * problem_obj.disc.D13(end - problem_obj.Lambda.negDimC + 1:end, :)]));
            
            
            obj.Boundary_Operators_a.H_W_Dy_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Aplus)  * problem_obj.disc.D21(1:problem_obj.Lambda.posDimA            , :); problem_obj.X.D21plus']));
            obj.Boundary_Operators_a.H_W_Dy_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Aminus) * problem_obj.disc.D21(end - problem_obj.Lambda.negDimA + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D21minus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim); problem_obj.X.D21minus']));

            obj.Boundary_Operators_a.H_E_Dy_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D21plus' ;-inv(problem_obj.Lambda.Aplus) * problem_obj.disc.D21(1:problem_obj.Lambda.posDimA            , :)]));
            obj.Boundary_Operators_a.H_E_Dy_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D21minus'; zeros(problem_obj.dim - size(problem_obj.X.D21minus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Aminus) * problem_obj.disc.D21(end - problem_obj.Lambda.negDimA + 1:end, :)]));

            obj.Boundary_Operators_a.H_S_Dy_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Bplus)  * problem_obj.disc.D22(1:problem_obj.Lambda.posDimB            , :); problem_obj.X.D22plus']));
            obj.Boundary_Operators_a.H_S_Dy_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Bminus) * problem_obj.disc.D22(end - problem_obj.Lambda.negDimB + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D22minus, 2) - length(problem_obj.Lambda.Bminus), problem_obj.dim); problem_obj.X.D22minus']));

            obj.Boundary_Operators_a.H_N_Dy_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D22plus' ;-inv(problem_obj.Lambda.Bplus) * problem_obj.disc.D22(1:problem_obj.Lambda.posDimB            , :)]));
            obj.Boundary_Operators_a.H_N_Dy_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D22minus'; zeros(problem_obj.dim - size(problem_obj.X.D22minus, 2) - length(problem_obj.Lambda.Bminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Bminus) * problem_obj.disc.D22(end - problem_obj.Lambda.negDimB + 1:end, :)]));

            obj.Boundary_Operators_a.H_B_Dy_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Cplus)  * problem_obj.disc.D23(1:problem_obj.Lambda.posDimC            , :); problem_obj.X.D23plus']));
            obj.Boundary_Operators_a.H_B_Dy_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Cminus) * problem_obj.disc.D23(end - problem_obj.Lambda.negDimC + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D23minus, 2) - length(problem_obj.Lambda.Cminus), problem_obj.dim); problem_obj.X.D23minus']));

            obj.Boundary_Operators_a.H_F_Dy_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D23plus' ;-inv(problem_obj.Lambda.Cplus) * problem_obj.disc.D23(1:problem_obj.Lambda.posDimC            , :)]));
            obj.Boundary_Operators_a.H_F_Dy_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D23minus'; zeros(problem_obj.dim - size(problem_obj.X.D23minus, 2) - length(problem_obj.Lambda.Cminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Cminus) * problem_obj.disc.D23(end - problem_obj.Lambda.negDimC + 1:end, :)]));

            
            obj.Boundary_Operators_a.H_W_Dz_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Aplus)  * problem_obj.disc.D31(1:problem_obj.Lambda.posDimA            , :); problem_obj.X.D31plus']));
            obj.Boundary_Operators_a.H_W_Dz_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Aminus) * problem_obj.disc.D31(end - problem_obj.Lambda.negDimA + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D31minus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim); problem_obj.X.D31minus']));

            obj.Boundary_Operators_a.H_E_Dz_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D31plus' ;-inv(problem_obj.Lambda.Aplus) * problem_obj.disc.D31(1:problem_obj.Lambda.posDimA            , :)]));
            obj.Boundary_Operators_a.H_E_Dz_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D31minus'; zeros(problem_obj.dim - size(problem_obj.X.D31minus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Aminus) * problem_obj.disc.D31(end - problem_obj.Lambda.negDimA + 1:end, :)]));

            obj.Boundary_Operators_a.H_S_Dz_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Bplus)  * problem_obj.disc.D32(1:problem_obj.Lambda.posDimB            , :); problem_obj.X.D32plus']));
            obj.Boundary_Operators_a.H_S_Dz_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Bminus) * problem_obj.disc.D32(end - problem_obj.Lambda.negDimB + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D32minus, 2) - length(problem_obj.Lambda.Bminus), problem_obj.dim); problem_obj.X.D32minus']));

            obj.Boundary_Operators_a.H_N_Dz_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D32plus' ;-inv(problem_obj.Lambda.Bplus) * problem_obj.disc.D32(1:problem_obj.Lambda.posDimB            , :)]));
            obj.Boundary_Operators_a.H_N_Dz_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D32minus'; zeros(problem_obj.dim - size(problem_obj.X.D32minus, 2) - length(problem_obj.Lambda.Bminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Bminus) * problem_obj.disc.D32(end - problem_obj.Lambda.negDimB + 1:end, :)]));

            obj.Boundary_Operators_a.H_B_Dz_p_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Cplus)  * problem_obj.disc.D33(1:problem_obj.Lambda.posDimC            , :); problem_obj.X.D33plus']));
            obj.Boundary_Operators_a.H_B_Dz_m_a = sparse(kron(Ixyz, em * [-inv(problem_obj.Lambda.Cminus) * problem_obj.disc.D33(end - problem_obj.Lambda.negDimC + 1:end, :); ...
                zeros(problem_obj.dim - size(problem_obj.X.D33minus, 2) - length(problem_obj.Lambda.Cminus), problem_obj.dim); problem_obj.X.D33minus']));

            obj.Boundary_Operators_a.H_F_Dz_p_a = sparse(kron(Ixyz, em * [problem_obj.X.D33plus' ;-inv(problem_obj.Lambda.Cplus) * problem_obj.disc.D33(1:problem_obj.Lambda.posDimC            , :)]));
            obj.Boundary_Operators_a.H_F_Dz_m_a = sparse(kron(Ixyz, em * [problem_obj.X.D33minus'; zeros(problem_obj.dim - size(problem_obj.X.D33minus, 2) - length(problem_obj.Lambda.Cminus), problem_obj.dim) ;...
                -inv(problem_obj.Lambda.Cminus) * problem_obj.disc.D33(end - problem_obj.Lambda.negDimC + 1:end, :)]));

            obj.Lambda.W_p_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.Aplus   , zeros(problem_obj.dim - length(problem_obj.Lambda.Aplus)  - length(problem_obj.Lambda.D11plus)),  problem_obj.Lambda.D11plus)));
            obj.Lambda.E_m_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.D11minus, zeros(problem_obj.dim - length(problem_obj.Lambda.Aminus) - length(problem_obj.Lambda.D11minus)), problem_obj.Lambda.Aminus)));
            obj.Lambda.S_p_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.Bplus   , zeros(problem_obj.dim - length(problem_obj.Lambda.Bplus)  - length(problem_obj.Lambda.D22plus)),  problem_obj.Lambda.D22plus)));
            obj.Lambda.N_m_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.D22minus, zeros(problem_obj.dim - length(problem_obj.Lambda.Bminus) - length(problem_obj.Lambda.D22minus)), problem_obj.Lambda.Bminus)));
            obj.Lambda.B_p_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.Cplus   , zeros(problem_obj.dim - length(problem_obj.Lambda.Cplus)  - length(problem_obj.Lambda.D33plus)),  problem_obj.Lambda.D33plus)));
            obj.Lambda.F_m_a = sparse(kron(Ixyz, blkdiag(problem_obj.Lambda.D33minus, zeros(problem_obj.dim - length(problem_obj.Lambda.Cminus) - length(problem_obj.Lambda.D33minus)), problem_obj.Lambda.Cminus)));
               
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
            Dxs   = kron(kron(kron(Dx, Iy), Iz), Isys);
            Dys   = kron(kron(kron(Ix, Dy), Iz), Isys);
            Dzs   = kron(kron(kron(Ix, Iy), Dz), Isys);

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
            
            H_E_Dx_p_n = sparse(Dxs * obj.Boundary_Operators_a.H_E_Dx_p_a);
            H_E_Dx_m_n = sparse(Dxs * obj.Boundary_Operators_a.H_E_Dx_m_a);
            H_W_Dx_p_n = sparse(Dxs * obj.Boundary_Operators_a.H_W_Dx_p_a);
            H_W_Dx_m_n = sparse(Dxs * obj.Boundary_Operators_a.H_W_Dx_m_a);
            H_N_Dx_p_n = sparse(Dxs * obj.Boundary_Operators_a.H_N_Dx_p_a);
            H_N_Dx_m_n = sparse(Dxs * obj.Boundary_Operators_a.H_N_Dx_m_a);
            H_S_Dx_p_n = sparse(Dxs * obj.Boundary_Operators_a.H_S_Dx_p_a);
            H_S_Dx_m_n = sparse(Dxs * obj.Boundary_Operators_a.H_S_Dx_m_a);
            H_F_Dx_p_n = sparse(Dxs * obj.Boundary_Operators_a.H_F_Dx_p_a);
            H_F_Dx_m_n = sparse(Dxs * obj.Boundary_Operators_a.H_F_Dx_m_a);
            H_B_Dx_p_n = sparse(Dxs * obj.Boundary_Operators_a.H_B_Dx_p_a);
            H_B_Dx_m_n = sparse(Dxs * obj.Boundary_Operators_a.H_B_Dx_m_a);
            
            H_E_Dy_p_n = sparse(Dys * obj.Boundary_Operators_a.H_E_Dy_p_a);
            H_E_Dy_m_n = sparse(Dys * obj.Boundary_Operators_a.H_E_Dy_m_a);
            H_W_Dy_p_n = sparse(Dys * obj.Boundary_Operators_a.H_W_Dy_p_a);
            H_W_Dy_m_n = sparse(Dys * obj.Boundary_Operators_a.H_W_Dy_m_a);
            H_N_Dy_p_n = sparse(Dys * obj.Boundary_Operators_a.H_N_Dy_p_a);
            H_N_Dy_m_n = sparse(Dys * obj.Boundary_Operators_a.H_N_Dy_m_a);
            H_S_Dy_p_n = sparse(Dys * obj.Boundary_Operators_a.H_S_Dy_p_a);
            H_S_Dy_m_n = sparse(Dys * obj.Boundary_Operators_a.H_S_Dy_m_a);
            H_F_Dy_p_n = sparse(Dys * obj.Boundary_Operators_a.H_F_Dy_p_a);
            H_F_Dy_m_n = sparse(Dys * obj.Boundary_Operators_a.H_F_Dy_m_a);
            H_B_Dy_p_n = sparse(Dys * obj.Boundary_Operators_a.H_B_Dy_p_a);
            H_B_Dy_m_n = sparse(Dys * obj.Boundary_Operators_a.H_B_Dy_m_a);
           
            H_E_Dz_p_n = sparse(Dzs * obj.Boundary_Operators_a.H_E_Dz_p_a);
            H_E_Dz_m_n = sparse(Dzs * obj.Boundary_Operators_a.H_E_Dz_m_a);
            H_W_Dz_p_n = sparse(Dzs * obj.Boundary_Operators_a.H_W_Dz_p_a);
            H_W_Dz_m_n = sparse(Dzs * obj.Boundary_Operators_a.H_W_Dz_m_a);
            H_N_Dz_p_n = sparse(Dzs * obj.Boundary_Operators_a.H_N_Dz_p_a);
            H_N_Dz_m_n = sparse(Dzs * obj.Boundary_Operators_a.H_N_Dz_m_a);
            H_S_Dz_p_n = sparse(Dzs * obj.Boundary_Operators_a.H_S_Dz_p_a);
            H_S_Dz_m_n = sparse(Dzs * obj.Boundary_Operators_a.H_S_Dz_m_a);
            H_F_Dz_p_n = sparse(Dzs * obj.Boundary_Operators_a.H_F_Dz_p_a);
            H_F_Dz_m_n = sparse(Dzs * obj.Boundary_Operators_a.H_F_Dz_m_a);
            H_B_Dz_p_n = sparse(Dzs * obj.Boundary_Operators_a.H_B_Dz_p_a);
            H_B_Dz_m_n = sparse(Dzs * obj.Boundary_Operators_a.H_B_Dz_m_a);
            
            obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n + H_E_Dx_p_n + H_E_Dy_p_n + H_E_Dz_p_n);
            obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n + H_E_Dx_m_n + H_E_Dy_m_n + H_E_Dz_m_n);
            obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n + H_W_Dx_p_n + H_W_Dy_p_n + H_W_Dz_p_n);
            obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n + H_W_Dx_m_n + H_W_Dy_m_n + H_W_Dz_m_n);
            obj.Boundary_Operators_n.H_N_p_n = sparse(H_N_0_p_n + H_N_Dx_p_n + H_N_Dy_p_n + H_N_Dz_p_n);
            obj.Boundary_Operators_n.H_N_m_n = sparse(H_N_0_m_n + H_N_Dx_m_n + H_N_Dy_m_n + H_N_Dz_m_n);
            obj.Boundary_Operators_n.H_S_p_n = sparse(H_S_0_p_n + H_S_Dx_p_n + H_S_Dy_p_n + H_S_Dz_p_n);
            obj.Boundary_Operators_n.H_S_m_n = sparse(H_S_0_m_n + H_S_Dx_m_n + H_S_Dy_m_n + H_S_Dz_m_n);
            obj.Boundary_Operators_n.H_F_p_n = sparse(H_F_0_p_n + H_F_Dx_p_n + H_F_Dy_p_n + H_F_Dz_p_n);
            obj.Boundary_Operators_n.H_F_m_n = sparse(H_F_0_m_n + H_F_Dx_m_n + H_F_Dy_m_n + H_F_Dz_m_n);
            obj.Boundary_Operators_n.H_B_p_n = sparse(H_B_0_p_n + H_B_Dx_p_n + H_B_Dy_p_n + H_B_Dz_p_n);
            obj.Boundary_Operators_n.H_B_m_n = sparse(H_B_0_m_n + H_B_Dx_m_n + H_B_Dy_m_n + H_B_Dz_m_n);

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
            
            X   = grid_obj.grid(:, :, :, 1);
            Y   = grid_obj.grid(:, :, :, 2);
            Z   = grid_obj.grid(:, :, :, 3);
            dim = problem_obj.dim;
            
            obj.Solution.u_a = reshape(reshape(problem_obj.data.u_a(X(:), Y(:), Z(:), t), [], dim)', [], 1);
                
        end
    end
end