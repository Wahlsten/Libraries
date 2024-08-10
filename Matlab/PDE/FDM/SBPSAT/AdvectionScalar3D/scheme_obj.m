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
        
        % Numerical solution
        Solution;               % Structure with solutions: u_a, u_n, u_ns, e_L2
        
    end
    methods
        function obj = createAMatrix(obj, grid_obj, sbp_obj, problem_obj)
            % Computes matrix A in : Ut + AU = Pen
            
            E_0            = cell(1, 3);
            E_N            = cell(1, 3);
            obj.Matrices.I = cell(1, 3);
            
            for k = 1:3 %grid_obj.dim
                
                E_0{k}            = sparse(grid_obj.N(k) + 1, grid_obj.N(k) + 1);
                E_N{k}            = sparse(grid_obj.N(k) + 1, grid_obj.N(k) + 1);
                E_0{k}(1, 1)      = 1;
                E_N{k}(end, end)  = 1;
                obj.Matrices.I{k} = speye(grid_obj.N(k) + 1);
                
            end
            
            obj.Scheme.E_East  = sparse(kron(kron(E_N{1}, obj.Matrices.I{2}), obj.Matrices.I{3}));
            obj.Scheme.E_West  = sparse(kron(kron(E_0{1}, obj.Matrices.I{2}), obj.Matrices.I{3}));
            obj.Scheme.E_North = sparse(kron(kron(obj.Matrices.I{1}, E_N{2}), obj.Matrices.I{3}));
            obj.Scheme.E_South = sparse(kron(kron(obj.Matrices.I{1}, E_0{2}), obj.Matrices.I{3}));
            obj.Scheme.E_Front = sparse(kron(kron(obj.Matrices.I{1}, obj.Matrices.I{2}), E_N{3}));
            obj.Scheme.E_Back  = sparse(kron(kron(obj.Matrices.I{1}, obj.Matrices.I{2}), E_0{3}));

            obj.Scheme.PinvENSigmaE = sparse(kron(kron(sbp_obj.Pinv{1}, obj.Matrices.I{2}), obj.Matrices.I{3}) * obj.Sigma.E * obj.Scheme.E_East);
            obj.Scheme.PinvE0SigmaW = sparse(kron(kron(sbp_obj.Pinv{1}, obj.Matrices.I{2}), obj.Matrices.I{3}) * obj.Sigma.W * obj.Scheme.E_West);
            obj.Scheme.PinvENSigmaN = sparse(kron(kron(obj.Matrices.I{1}, sbp_obj.Pinv{2}), obj.Matrices.I{3}) * obj.Sigma.N * obj.Scheme.E_North);
            obj.Scheme.PinvE0SigmaS = sparse(kron(kron(obj.Matrices.I{1}, sbp_obj.Pinv{2}), obj.Matrices.I{3}) * obj.Sigma.S * obj.Scheme.E_South);
            obj.Scheme.PinvENSigmaF = sparse(kron(kron(obj.Matrices.I{1}, obj.Matrices.I{2}), sbp_obj.Pinv{3}) * obj.Sigma.F * obj.Scheme.E_Front);
            obj.Scheme.PinvE0SigmaB = sparse(kron(kron(obj.Matrices.I{1}, obj.Matrices.I{2}), sbp_obj.Pinv{3}) * obj.Sigma.B * obj.Scheme.E_Back);

            obj.Scheme.PinvE0SigmaWH0 = sparse(obj.Scheme.PinvE0SigmaW * sparse(obj.Boundary_Operators_a.H_W_0_p_a));
            obj.Scheme.PinvENSigmaEH0 = sparse(obj.Scheme.PinvENSigmaE * sparse(obj.Boundary_Operators_a.H_E_0_m_a));
            obj.Scheme.PinvE0SigmaSH0 = sparse(obj.Scheme.PinvE0SigmaS * sparse(obj.Boundary_Operators_a.H_S_0_p_a));
            obj.Scheme.PinvENSigmaNH0 = sparse(obj.Scheme.PinvENSigmaN * sparse(obj.Boundary_Operators_a.H_N_0_m_a));
            obj.Scheme.PinvE0SigmaBH0 = sparse(obj.Scheme.PinvE0SigmaB * sparse(obj.Boundary_Operators_a.H_B_0_p_a));
            obj.Scheme.PinvENSigmaFH0 = sparse(obj.Scheme.PinvENSigmaF * sparse(obj.Boundary_Operators_a.H_F_0_m_a));
            
            obj.Scheme.Amat    = problem_obj.cont.a;
            obj.Scheme.Bmat    = problem_obj.cont.b;
            obj.Scheme.Cmat    = problem_obj.cont.c;

            Am  = obj.Scheme.Amat;
            Bm  = obj.Scheme.Bmat;
            Cm  = obj.Scheme.Cmat;

            Dx = sparse(kron(kron(sbp_obj.D1{1}    , obj.Matrices.I{2}), obj.Matrices.I{3}));
            Dy = sparse(kron(kron(obj.Matrices.I{1}, sbp_obj.D1{2})    , obj.Matrices.I{3}));
            Dz = sparse(kron(kron(obj.Matrices.I{1}, obj.Matrices.I{2}), sbp_obj.D1{3}));

            obj.A = - (Am * Dx + Bm * Dy + Cm * Dz) ...
                + obj.Scheme.PinvE0SigmaW * obj.Boundary_Operators_n.W ...
                + obj.Scheme.PinvENSigmaE * obj.Boundary_Operators_n.E ...
                + obj.Scheme.PinvE0SigmaS * obj.Boundary_Operators_n.S ...
                + obj.Scheme.PinvENSigmaN * obj.Boundary_Operators_n.N ...
                + obj.Scheme.PinvENSigmaF * obj.Boundary_Operators_n.F ...
                + obj.Scheme.PinvE0SigmaB * obj.Boundary_Operators_n.B;

            obj.A = sparse(obj.A);

        end
        function obj = createScheme(obj, grid_obj, problem_obj)
            
            X   = grid_obj.grid(:, :, :, 1);
            Y   = grid_obj.grid(:, :, :, 2);
            Z   = grid_obj.grid(:, :, :, 3);
            
            A_a = problem_obj.cont.a;
            B_a = problem_obj.cont.b;
            C_a = problem_obj.cont.c;

            force = @(t) (   reshape(problem_obj.data.uT_a (X(:), Y(:), Z(:), t)', [], 1) ...
                +  A_a     * reshape(problem_obj.data.uX_a (X(:), Y(:), Z(:), t)', [], 1) ...
                +  B_a     * reshape(problem_obj.data.uY_a (X(:), Y(:), Z(:), t)', [], 1) ...
                +  C_a     * reshape(problem_obj.data.uZ_a (X(:), Y(:), Z(:), t)', [], 1));

            obj.cA = @(t, u) (obj.A * u ...
                - obj.Scheme.PinvE0SigmaWH0  * reshape(problem_obj.data.gW  (X(:), Y(:), Z(:), t)', [], 1) ...
                - obj.Scheme.PinvENSigmaEH0  * reshape(problem_obj.data.gE  (X(:), Y(:), Z(:), t)', [], 1) ...
                - obj.Scheme.PinvE0SigmaSH0  * reshape(problem_obj.data.gS  (X(:), Y(:), Z(:), t)', [], 1) ...
                - obj.Scheme.PinvENSigmaNH0  * reshape(problem_obj.data.gN  (X(:), Y(:), Z(:), t)', [], 1) ...
                - obj.Scheme.PinvENSigmaFH0  * reshape(problem_obj.data.gF  (X(:), Y(:), Z(:), t)', [], 1) ...
                - obj.Scheme.PinvE0SigmaBH0  * reshape(problem_obj.data.gB  (X(:), Y(:), Z(:), t)', [], 1) ...
                + force(t));

        end
        function obj = createPenalties(obj)

            obj.Sigma.E =  obj.Boundary_Operators_n.H_E_m_n' * obj.Lambda.E_m;
            obj.Sigma.W = -obj.Boundary_Operators_n.H_W_p_n' * obj.Lambda.W_p;
            obj.Sigma.N =  obj.Boundary_Operators_n.H_N_m_n' * obj.Lambda.N_m;
            obj.Sigma.S = -obj.Boundary_Operators_n.H_S_p_n' * obj.Lambda.S_p;
            obj.Sigma.F = -obj.Boundary_Operators_n.H_F_m_n' * obj.Lambda.F_m;
            obj.Sigma.B = -obj.Boundary_Operators_n.H_B_p_n' * obj.Lambda.B_p;
            
        end
        function obj = createContinuousBoundaryOperators(obj, problem_obj, grid_obj)
            
                
                Am = diag(problem_obj.cont.a);
                Bm = diag(problem_obj.cont.b);
                Cm = diag(problem_obj.cont.c);
                
                obj.Boundary_Operators_a.H_E_0_p_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am >  1e-12) .* Am));
                obj.Boundary_Operators_a.H_E_0_m_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am < -1e-12) .* Am));
                obj.Boundary_Operators_a.H_W_0_p_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am >  1e-12) .* Am));
                obj.Boundary_Operators_a.H_W_0_m_a  = sparse(diag( (abs(Am) < 1e-12) .* abs(Am) + (Am < -1e-12) .* Am));
                
                obj.Boundary_Operators_a.H_N_0_p_a  = sparse(diag( (abs(Bm) < 1e-12) .* abs(Bm) + (Bm >  1e-12) .* Bm));
                obj.Boundary_Operators_a.H_N_0_m_a  = sparse(diag( (abs(Bm) < 1e-12) .* abs(Bm) + (Bm < -1e-12) .* Bm));
                obj.Boundary_Operators_a.H_S_0_p_a  = sparse(diag( (abs(Bm) < 1e-12) .* abs(Bm) + (Bm >  1e-12) .* Bm));
                obj.Boundary_Operators_a.H_S_0_m_a  = sparse(diag( (abs(Bm) < 1e-12) .* abs(Bm) + (Bm < -1e-12) .* Bm));

                obj.Boundary_Operators_a.H_F_0_p_a  = sparse(diag( (abs(Cm) < 1e-12) .* abs(Cm) + (Cm >  1e-12) .* Cm));
                obj.Boundary_Operators_a.H_F_0_m_a  = sparse(diag( (abs(Cm) < 1e-12) .* abs(Cm) + (Cm < -1e-12) .* Cm));
                obj.Boundary_Operators_a.H_B_0_p_a  = sparse(diag( (abs(Cm) < 1e-12) .* abs(Cm) + (Cm >  1e-12) .* Cm));
                obj.Boundary_Operators_a.H_B_0_m_a  = sparse(diag( (abs(Cm) < 1e-12) .* abs(Cm) + (Cm < -1e-12) .* Cm));

                obj.Lambda.E_m_a = sparse(diag(-pinv(diag(abs(full(Am)))) * (abs(Am) > 1e-12)));
                obj.Lambda.W_p_a = sparse(diag( pinv(diag(abs(full(Am)))) * (abs(Am) > 1e-12)));
                obj.Lambda.N_m_a = sparse(diag(-pinv(diag(abs(full(Bm)))) * (abs(Bm) > 1e-12)));
                obj.Lambda.S_p_a = sparse(diag( pinv(diag(abs(full(Bm)))) * (abs(Bm) > 1e-12)));                    
                obj.Lambda.F_m_a = sparse(diag(-pinv(diag(abs(full(Cm)))) * (abs(Cm) > 1e-12)));
                obj.Lambda.B_p_a = sparse(diag( pinv(diag(abs(full(Cm)))) * (abs(Cm) > 1e-12)));

        end
        function obj = createBoundaryOperator(obj, problem_obj, grid_obj)
            
                obj = createContinuousBoundaryOperators(obj, problem_obj, grid_obj);

                H_E_0_p_n = sparse(obj.Boundary_Operators_a.H_E_0_p_a);
                H_E_0_m_n = sparse(obj.Boundary_Operators_a.H_E_0_m_a);
                H_W_0_p_n = sparse(obj.Boundary_Operators_a.H_W_0_p_a);
                H_W_0_m_n = sparse(obj.Boundary_Operators_a.H_W_0_m_a);
                H_N_0_p_n = sparse(obj.Boundary_Operators_a.H_N_0_p_a);
                H_N_0_m_n = sparse(obj.Boundary_Operators_a.H_N_0_m_a);
                H_S_0_p_n = sparse(obj.Boundary_Operators_a.H_S_0_p_a);
                H_S_0_m_n = sparse(obj.Boundary_Operators_a.H_S_0_m_a);
                H_F_0_p_n = sparse(obj.Boundary_Operators_a.H_F_0_p_a);
                H_F_0_m_n = sparse(obj.Boundary_Operators_a.H_F_0_m_a);
                H_B_0_p_n = sparse(obj.Boundary_Operators_a.H_B_0_p_a);
                H_B_0_m_n = sparse(obj.Boundary_Operators_a.H_B_0_m_a);
                
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
                
                obj.Lambda.E_m = sparse(diag(-pinv(abs(full(diag(problem_obj.cont.a)))) ...
                    * (abs(diag(problem_obj.cont.a)) > 1e-12)));
                obj.Lambda.W_p = sparse(diag( pinv(abs(full(diag(problem_obj.cont.a)))) ...
                    * (abs(diag(problem_obj.cont.a)) > 1e-12)));
                obj.Lambda.N_m = sparse(diag(-pinv(abs(full(diag(problem_obj.cont.b)))) ...
                    * (abs(diag(problem_obj.cont.b)) > 1e-12)));
                obj.Lambda.S_p = sparse(diag( pinv(abs(full(diag(problem_obj.cont.b)))) ...
                    * (abs(diag(problem_obj.cont.b)) > 1e-12)));
                obj.Lambda.F_m = sparse(diag(-pinv(abs(full(diag(problem_obj.cont.c)))) ...
                    * (abs(diag(problem_obj.cont.c)) > 1e-12)));
                obj.Lambda.B_p = sparse(diag( pinv(abs(full(diag(problem_obj.cont.c)))) ...
                    * (abs(diag(problem_obj.cont.c)) > 1e-12)));

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

            opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

            [~, obj.Solution.u_n] = ode45(obj.cA, grid_obj.t, ...
                reshape(problem_obj.data.f(X(:), Y(:), Z(:), 0)', [], 1), opts);

            obj.Solution.u_n = obj.Solution.u_n';

        end
        function plotSolution(obj, grid_obj)
            
            figure
            plot(grid_obj.grid, obj.Solution.u_n(:, 1))

        end
        function obj = computeError(obj, grid_obj, sbp_obj, problem_obj)
            
            obj   = computeAnalyticalSolution(obj, grid_obj, problem_obj, grid_obj.maxv(end));
            error = abs(obj.Solution.u_n(:,end) - obj.Solution.u_a);
            obj.Solution.e_L2 = sqrt(error' * kron(kron(inv(sbp_obj.Pinv{1}), inv(sbp_obj.Pinv{2})), inv(sbp_obj.Pinv{3})) * error);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, t)
            
            X = grid_obj.grid(:,:,:,1);
            Y = grid_obj.grid(:,:,:,2);
            Z = grid_obj.grid(:,:,:,3);
                
            obj.Solution.u_a = reshape(problem_obj.data.u_a(X(:), Y(:), Z(:), t)', [], 1);

        end
    end
end