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
            
            E_0               = zeros(grid_obj.N(1) + 1);
            E_N               = zeros(grid_obj.N(1) + 1);
            E_0(1, 1)         = 1;
            E_N(end, end)     = 1;
            obj.Matrices.I    = eye(grid_obj.N(1) + 1);
            obj.Matrices.Isys = eye(problem_obj.dim);
            
            obj.Scheme.E_East = sparse(kron(E_N, obj.Matrices.Isys));
            obj.Scheme.E_West = sparse(kron(E_0, obj.Matrices.Isys));
            
            obj.Scheme.PinvENSigmaE = sparse(kron(sbp_obj.Pinv, obj.Matrices.Isys) * obj.Sigma.E * obj.Scheme.E_East);
            obj.Scheme.PinvE0SigmaW = sparse(kron(sbp_obj.Pinv, obj.Matrices.Isys) * obj.Sigma.W * obj.Scheme.E_West);
            
            obj.Scheme.PinvE0SigmaWH0  = sparse(obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_a.H_W_0_p_a);
            obj.Scheme.PinvE0SigmaWHDx = sparse(obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_a.H_W_Dx_p_a);
            obj.Scheme.PinvENSigmaEH0  = sparse(obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_a.H_E_0_m_a);
            obj.Scheme.PinvENSigmaEHDx = sparse(obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_a.H_E_Dx_m_a);
            
            obj.Scheme.Amat    = obj.Trans.A;
            obj.Scheme.Bmat    = obj.Trans.B;
            obj.Scheme.epsilon = problem_obj.epsilon;

            Am  = obj.Scheme.Amat;
            Bm  = obj.Scheme.Bmat;
            em  = obj.Scheme.epsilon;

            obj.Scheme.epsBmat = sparse(problem_obj.epsilon * problem_obj.disc.B);

            Dx = sparse(kron(sbp_obj.D1, obj.Matrices.Isys));

            obj.A = - (Am * Dx) ...
                +      em * Bm * Dx * Dx          ...
                + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E;

            obj.A = sparse(obj.A);

        end
        function obj = createScheme(obj, grid_obj, problem_obj)
            
            X = grid_obj.grid;
                            
                xlen     = grid_obj.N(1) + 1;
                xzeroVec = grid_obj.minv(1) * ones(1, xlen);
                xoneVec  = grid_obj.maxv(1) * ones(1, xlen);
                
                if strcmp(problem_obj.type.coeff, 'constant')

                    Amat     = sparse(kron(obj.Matrices.I, problem_obj.disc.A));
                    epsBmat  = sparse(kron(obj.Matrices.I, problem_obj.epsilon * problem_obj.disc.B));

                elseif strcmp(problem_obj.type.coeff, 'variable')

                    Amat     = sparse(problem_obj.disc.A);
                    epsBmat  = sparse(problem_obj.epsilon * problem_obj.disc.B);

                end

                force = @(t) (   reshape(problem_obj.data.uT_a (X, t), [], 1) ...
                    +  Amat    * reshape(problem_obj.data.uX_a (X, t), [], 1) ...
                    -  epsBmat * reshape(problem_obj.data.uXX_a(X, t), [], 1));

                obj.cA = @(t, u) (obj.A * u ...
                    - obj.Scheme.PinvE0SigmaWH0  * reshape(problem_obj.data.gW (xzeroVec, t), [], 1) ...
                    - obj.Scheme.PinvE0SigmaWHDx * reshape(problem_obj.data.gWD(xzeroVec, t), [], 1) ...
                    - obj.Scheme.PinvENSigmaEH0  * reshape(problem_obj.data.gE (xoneVec,  t), [], 1) ...
                    - obj.Scheme.PinvENSigmaEHDx * reshape(problem_obj.data.gED(xoneVec,  t), [], 1) ...
                    + force(t));
                        
        end
        function obj = createPenalties(obj)

            obj.Sigma.E =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            
        end
        function obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj)
            
            if strcmp(problem_obj.type.coeff, 'constant')
                
                Ix   = eye(grid_obj.N(1) + 1);
                
                obj.Boundary_Operators_a.H_W_0_p_a  = sparse(kron(Ix, blkdiag(eye(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), zeros(problem_obj.Lambda.negDimA))));
                obj.Boundary_Operators_a.H_W_0_m_a  = sparse(kron(Ix, diag(ones(1, problem_obj.Lambda.negDimA), problem_obj.Lambda.posDimA + problem_obj.Lambda.zeroDimA)));
                
                obj.Boundary_Operators_a.H_W_Dx_p_a = sparse(kron(Ix, problem_obj.epsilon * [-inv(problem_obj.Lambda.Aplus)  * problem_obj.disc.B(1:problem_obj.Lambda.posDimA            , :); problem_obj.X.Bplus']));
                obj.Boundary_Operators_a.H_W_Dx_m_a = sparse(kron(Ix, problem_obj.epsilon * [-inv(problem_obj.Lambda.Aminus) * problem_obj.disc.B(end - problem_obj.Lambda.negDimA + 1:end, :); ...
                    zeros(problem_obj.dim - size(problem_obj.X.Bminus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim); problem_obj.X.Bminus']));
                
                obj.Boundary_Operators_a.H_E_0_p_a  = sparse(kron(Ix, diag([ones(1, problem_obj.Lambda.posDimA), zeros(1, problem_obj.Lambda.zeroDimA), zeros(1, problem_obj.Lambda.negDimA)])'));
                obj.Boundary_Operators_a.H_E_0_m_a  = sparse(kron(Ix, blkdiag(zeros(problem_obj.Lambda.posDimA), zeros(problem_obj.Lambda.zeroDimA), eye(problem_obj.Lambda.negDimA))));
                
                obj.Boundary_Operators_a.H_E_Dx_p_a = sparse(kron(Ix, problem_obj.epsilon * [problem_obj.X.Bplus' ;-inv(problem_obj.Lambda.Aplus) * problem_obj.disc.B(1:problem_obj.Lambda.posDimA            , :)]));
                obj.Boundary_Operators_a.H_E_Dx_m_a = sparse(kron(Ix, problem_obj.epsilon * [problem_obj.X.Bminus'; zeros(problem_obj.dim - size(problem_obj.X.Bminus, 2) - length(problem_obj.Lambda.Aminus), problem_obj.dim) ;...
                    -inv(problem_obj.Lambda.Aminus) * problem_obj.disc.B(end - problem_obj.Lambda.negDimA + 1:end, :)]));
                
                obj.Lambda.W_p_a = sparse(kron(Ix, blkdiag(problem_obj.Lambda.Aplus , zeros(problem_obj.dim - length(problem_obj.Lambda.Aplus)  - length(problem_obj.Lambda.Bplus)),  problem_obj.Lambda.Bplus)));
                obj.Lambda.E_m_a = sparse(kron(Ix, blkdiag(problem_obj.Lambda.Bminus, zeros(problem_obj.dim - length(problem_obj.Lambda.Aminus) - length(problem_obj.Lambda.Bminus)), problem_obj.Lambda.Aminus)));
                
            elseif strcmp(problem_obj.type.coeff, 'variable')
                
                Am = diag(obj.Trans.A);
                Bm = diag(obj.Trans.B);
                
                obj.Boundary_Operators_a.H_E_0_p_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am >  1e-12) .* Am));
                obj.Boundary_Operators_a.H_E_0_m_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am < -1e-12) .* Am));
                obj.Boundary_Operators_a.H_W_0_p_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am >  1e-12) .* Am));
                obj.Boundary_Operators_a.H_W_0_m_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am < -1e-12) .* Am));
                
                obj.Boundary_Operators_a.H_E_Dx_p_a  = sparse(diag( (abs(Am) < 1e-12) .* -problem_obj.epsilon + -problem_obj.epsilon .* (abs(Am) > 1e-12) .* Bm));
                obj.Boundary_Operators_a.H_E_Dx_m_a  = sparse(diag( (abs(Am) < 1e-12) .* -problem_obj.epsilon + -problem_obj.epsilon .* (abs(Am) > 1e-12) .* Bm));
                obj.Boundary_Operators_a.H_W_Dx_p_a  = sparse(diag(-(abs(Am) < 1e-12) .*  problem_obj.epsilon + -problem_obj.epsilon .* (abs(Am) > 1e-12) .* Bm));
                obj.Boundary_Operators_a.H_W_Dx_m_a  = sparse(diag(-(abs(Am) < 1e-12) .*  problem_obj.epsilon + -problem_obj.epsilon .* (abs(Am) > 1e-12) .* Bm));
                
                obj.Lambda.E_m_a = sparse(diag(-(abs(Am) < 1e-12) .* problem_obj.epsilon .* Bm + -pinv(diag(abs(full(Am)))) * (abs(Am) > 1e-12)));
                obj.Lambda.W_p_a = sparse(diag( (abs(Am) < 1e-12) .* problem_obj.epsilon .* Bm +  pinv(diag(abs(full(Am)))) * (abs(Am) > 1e-12)));
                
            end
        end
        function obj = createBoundaryOperator(obj, grid_obj, sbp_obj, problem_obj)

            obj = transformation(obj, grid_obj, problem_obj);
            obj = createContinuousBoundaryOperators(obj, grid_obj, problem_obj);

            if strcmp(problem_obj.type.coeff, 'variable')

                Dx = sbp_obj.D1{1};
                Dy = sbp_obj.D1{2};
                Ix = eye(size(Dx));
                Iy = eye(size(Dy));

                H_E_0_p_n = sparse(obj.Boundary_Operators_a.H_E_0_p_a);
                H_E_0_m_n = sparse(obj.Boundary_Operators_a.H_E_0_m_a);
                H_W_0_p_n = sparse(obj.Boundary_Operators_a.H_W_0_p_a);
                H_W_0_m_n = sparse(obj.Boundary_Operators_a.H_W_0_m_a);
                
                H_E_Dx_p_n = sparse(obj.Boundary_Operators_a.H_E_Dx_p_a * kron(Dx, Iy));
                H_E_Dx_m_n = sparse(obj.Boundary_Operators_a.H_E_Dx_m_a * kron(Dx, Iy));
                H_W_Dx_p_n = sparse(obj.Boundary_Operators_a.H_W_Dx_p_a * kron(Dx, Iy));
                H_W_Dx_m_n = sparse(obj.Boundary_Operators_a.H_W_Dx_m_a * kron(Dx, Iy));
                
                obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n + H_E_Dx_p_n + H_E_Dy_p_n);
                obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n + H_E_Dx_m_n + H_E_Dy_m_n);
                obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n + H_W_Dx_p_n + H_W_Dy_p_n);
                obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n + H_W_Dx_m_n + H_W_Dy_m_n);
                
                obj.Lambda.E_m = sparse(diag(-pinv(abs(full(diag(obj.Trans.A)))) ...
                    * (abs(diag(obj.Trans.A)) > 1e-12)));
                obj.Lambda.W_p = sparse(diag( pinv(abs(full(diag(obj.Trans.A)))) ...
                    * (abs(diag(obj.Trans.A)) > 1e-12)));
                
            elseif strcmp(problem_obj.type.coeff, 'constant')

                Isys = eye(problem_obj.dim);
                Dx   = sbp_obj.D1;
                Ix   = eye(size(Dx));

                H_E_0_p_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_E_0_p_a);
                H_E_0_m_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_E_0_m_a);
                H_W_0_p_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_W_0_p_a);
                H_W_0_m_n = sparse(kron(Ix, Isys) * obj.Boundary_Operators_a.H_W_0_m_a);
                
                H_E_Dx_p_n = sparse(kron(Dx, Isys) * obj.Boundary_Operators_a.H_E_Dx_p_a);
                H_E_Dx_m_n = sparse(kron(Dx, Isys) * obj.Boundary_Operators_a.H_E_Dx_m_a);
                H_W_Dx_p_n = sparse(kron(Dx, Isys) * obj.Boundary_Operators_a.H_W_Dx_p_a);
                H_W_Dx_m_n = sparse(kron(Dx, Isys) * obj.Boundary_Operators_a.H_W_Dx_m_a);
                
                obj.Boundary_Operators_n.H_E_p_n = sparse(H_E_0_p_n + H_E_Dx_p_n);
                obj.Boundary_Operators_n.H_E_m_n = sparse(H_E_0_m_n + H_E_Dx_m_n);
                obj.Boundary_Operators_n.H_W_p_n = sparse(H_W_0_p_n + H_W_Dx_p_n);
                obj.Boundary_Operators_n.H_W_m_n = sparse(H_W_0_m_n + H_W_Dx_m_n);
                
                obj.Lambda.E_m = sparse(obj.Lambda.E_m_a);
                obj.Lambda.W_p = sparse(obj.Lambda.W_p_a);
                
            end

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
                                                
            if strcmp(problem_obj.type.coeff, 'variable')

                obj.Trans.A    = sparse(problem_obj.disc.A);
                obj.Trans.B    = sparse(problem_obj.disc.B);

            elseif strcmp(problem_obj.type.coeff, 'constant')

                Ix = eye(grid_obj.N(1) + 1, grid_obj.N(1) + 1);

                obj.Trans.A    = sparse(kron(Ix, problem_obj.disc.A));
                obj.Trans.B    = sparse(kron(Ix, problem_obj.disc.B));

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
            obj.Solution.e_L2 = sqrt(error' * kron(inv(sbp_obj.Pinv), obj.Matrices.Isys) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = grid_obj.grid;
            obj.Solution.u_a = reshape(problem_obj.data.u_a(X, t), [], 1);
            
        end
        
    end
end