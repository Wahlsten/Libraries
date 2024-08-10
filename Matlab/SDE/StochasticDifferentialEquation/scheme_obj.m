classdef scheme_obj
    % Numerical formulation
    % ut + Au = 0
    properties        
        % Numerical solution
        Solution;               % Structure with solutions: u_a, u_n, u_ns, e_L2
        
    end
    methods
        function obj = solveUsingEuler(obj, grid_obj, stc_obj, problem_obj)
            
            Z      = stc_obj.grid;
            a      = problem_obj.cont.a;
            b      = problem_obj.cont.b;
            b_x    = problem_obj.cont.b_x;
            R      = stc_obj.R;
            dt     = stc_obj.R * grid_obj.d * ones(stc_obj.N_mean, 1);
            X      = zeros(stc_obj.N_mean, grid_obj.N/R);
            X(:,1) = problem_obj.data.f(0);
            
            
            if strcmp(stc_obj.method, 'Euler')
            
                for k = 1:grid_obj.N/R
                    
                    tp = stc_obj.grid(R*k);
                    Xp = X(:, k);
                    Zp = cumsum(Z(:, R*(k-1) + 1:R*k), 2);
                    X(:, k + 1) = X(:, k) + a(Xp, tp) .* dt + b(Xp, tp) .* Zp(:, end);

                end

            elseif strcmp(stc_obj.method, 'Milstein')

                for k = 1:grid_obj.N/R

                    tp = stc_obj.grid(R*k);
                    Xp = X(:, k);
                    Zp = cumsum(Z(:, R*(k-1) + 1:R*k), 2);
                    X(:, k + 1) = X(:, k) + a(Xp, tp) .* dt + b(Xp, tp) .* Zp(:, end) ...
                        + 1/2 * b_x(Xp, tp) .* (Zp(:, end).^2 - dt);

                end
                
            elseif strcmp(stc_obj.method, 'RK')
                
                for k = 1:grid_obj.N/R
                    
                    tp = stc_obj.grid(R*k);
                    Xp = X(:, k);
                    Zp = cumsum(Z(:, R*(k-1) + 1:R*k), 2);
                    Yp = Xp + a(Xp, tp) .* dt + b(Xp, tp) .* sqrt(dt);
                    
                    X(:, k + 1) = X(:, k) + a(Xp, tp) .* dt + b(Xp, tp) .* Zp(:, end) ...
                        + 1/2 * (b(Yp, tp) - b(Xp, tp)) .* (Zp(:, end).^2 - dt) * 1./sqrt(dt);
                    
                end
                
            elseif strcmp(stc_obj.method, 'SBPinTime')
                
                acc     = 2;
                sbp     = sbp_obj;
                sbp.acc = acc;
                sbp     = createOperators(sbp, grid_obj);
                
                Sigma              = -1;
                Scheme.E_0         = sparse(zeros(grid_obj.N(1) + 1));
                Scheme.E_0(1,1)    = 1;
                Scheme.PinvE0Sigma = sparse(sbp.Pinv * Sigma * Scheme.E_0);
                f                  = problem_obj.data.f(grid_obj.t');
                Dt                 = sparse(sbp.D1);
                
                force = a(0, grid_obj.t') + b(0, grid_obj.t') .* [0; Z'];

                A = sparse(Dt - Scheme.PinvE0Sigma);
                F = force - Scheme.PinvE0Sigma * f;
                
                X = A\F;
                
            end
            
            obj.Solution.u_n = X;
            
        end
        function plotSolution(obj, grid_obj, stc_obj)

            R      = stc_obj.R;
            t      = grid_obj.t(1:R:end);
            
            figure
            hold on
            plot(t, mean(obj.Solution.u_n, 1))
            plot(t, mean(obj.Solution.u_a, 1))
            hold off
            
            legend('u_n', 'u_a')
            
%             figure
%             plot(grid_obj.t, obj.Solution.error)
%             
        end
        function obj = computeError(obj, grid_obj, problem_obj, stc_obj)
            
            obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, stc_obj);

            obj.Solution.error_strong = mean(abs(obj.Solution.u_a(:,end)  - obj.Solution.u_n(:,end)), 1);
            obj.Solution.error_weak   = abs(mean(obj.Solution.u_n(:,end)) - obj.Solution.mean_u_a);
            
        end
        function obj = computeAnalyticalSolution(obj, grid_obj, problem_obj, stc_obj)
            
            X0     = problem_obj.data.f(0) * ones(stc_obj.N_mean, 1);
            t      = grid_obj.t;
            dW     = stc_obj.grid;
            W      = cumsum(dW, 2);
            a      = problem_obj.cont.a(X0, t(2:end));
            b      = problem_obj.cont.b(X0, t(2:end));
            
            obj.Solution.u_a = X0 .* exp((a - 0.5 * b.^2) .* kron(ones(stc_obj.N_mean, 1), t(2:end)) + b .* W);
            
            obj.Solution.u_a = [X0, obj.Solution.u_a];
            
            obj.Solution.mean_u_a = X0(end) * exp(problem_obj.cont.a(X0(end), t(end)));
            
            % obj.Solution.u_a = X0 .* exp((a - 0.5 * b.^2) .* kron(ones(stc_obj.N_mean, 1), t(2:end)) ...
            %+ b .* W(:, 1:R:end));

        end
    end
end