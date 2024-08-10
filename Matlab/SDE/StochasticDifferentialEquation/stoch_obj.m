classdef stoch_obj
    % Uncertainty quantification object
    properties
        Distribution;
        grid;
        method;
        N_mean;
        dw;
        R;
    end
    methods
        function obj      = createGridS(obj, grid_obj)
            randn('state',100)
            
%             for k = 1:obj.N_mean
%                 
%                 obj.grid(k,:) = randn(1, grid_obj.N);
%             
%             end
%             
%             obj.grid = sqrt(grid_obj.d) * obj.grid;

              obj.grid = sqrt(grid_obj.d) * randn(obj.N_mean, grid_obj.N);
%             obj.dw   = obj.grid(:, 2:end) - obj.grid(:, 1:end - 1);
            
            %obj.grid = [zeros(obj.N_mean, 1) , obj.grid];
            
        end
        function moment   = computeMoment(obj, u_ns, m)
            
        end
        function variance = computeVariance(obj, u_ns)
            
        end
        function obj      = computeStatistics(obj, u_ns)
            
        end
        function obj      = computeAnalyticalStatistics(obj, problem_obj)
            
        end
    end
end