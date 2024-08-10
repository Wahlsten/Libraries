classdef error_obj
    % Error object
    properties
        
        N;       % Matrix of number of grid points
        M;       % Number of realizations used as Average 
        R;       % Resolution numerical method
        method;  % Vector of order of accuracy
        t;       % Matrix of times
        Error;   % Matrix of L_2 errors
        Order;   % Matrix of computed order of accuracy
        NoExact; % No exact solution
        Acc;     % Accuracy of the method
        
    end
    methods
        function obj = initialize(obj, problem, data)
            
            obj.Error.strong = zeros(length(obj.method), size(obj.R, 2) - obj.NoExact);
            obj.Error.weak   = zeros(length(obj.method), size(obj.R, 2) - obj.NoExact);
            obj.Order.strong = zeros(length(obj.method), size(obj.R, 2) - 1 - obj.NoExact);
            obj.Order.weak   = zeros(length(obj.method), size(obj.R, 2) - 1 - obj.NoExact);
            obj.t     = [];
            format long
            
        end
        function Table(obj, str)
            
            xVec = obj.N(1, :);
            order = obj.Acc;
            
            if strcmp(str, 'Order strong')
                
                titleH    = ' Order Of Accuracy ';
                spacing   = '\t';
                xVec(end) = [];
                data      = obj.Order.strong;
                
            elseif strcmp(str, 'Order weak')
                
                titleH    = ' Order Of Accuracy ';
                spacing   = '\t';
                xVec(end) = [];
                data      = obj.Order.weak;
                
            elseif strcmp(str, 'Error strong')
                
                titleH  = ' Error ';
                spacing = '\t\t';
                data    = obj.Error.strong;
                
            elseif strcmp(str, 'Error weak')
                
                titleH  = ' Error ';
                spacing = '\t\t';
                data    = obj.Error.weak;
                
            elseif strcmp(str, 'Time')
                
                titleH  = ' Time (s)';
                spacing = '\t\t';
                data    = obj.t;
                
            end
            
            oneEl = 7.5;
            numEl = round(length(spacing) / 2 * oneEl*(length(xVec)) + 2*oneEl);
            procSigns = repmat('%',1,round((numEl - length(titleH))/2));
            titleH = [procSigns, titleH, procSigns];
            fprintf('\n')
            fprintf('\n')
            disp(titleH)
            fprintf('\n')
            fprintf('Order/Nx')
            fprintf('\t')
            
            for i = 1:length(xVec)
                
                fprintf(num2str(xVec(i)))
                fprintf(spacing)
                
            end
            
            fprintf('\n')
            fprintf(num2str(repmat('-', 1, numEl)))
            fprintf('\n')
            
            for k = 1:size(data, 1)
                
                fprintf(num2str(order(k)))
                fprintf('\t\t')
                
                for l = 1:size(data,2)
                    
                    fprintf('%12.8f', data(k,l)); %num2str(data(k,l)))
                    fprintf('\t')
                    
                end
                
                fprintf('\n')
                
            end
            
        end
        
        function obj = StrongConvergence(obj, problem, data)
            
            figure
            
            for m = 1:length(obj.method)
                            
                if obj.NoExact
                   
                    [s, ~] = SDESolver(problem, data, obj.N(end), 0, 1, obj.method{m}, obj.R(1), obj.M);
                    
                    obj.Error.u_n(m, 1) = mean(s.Solution.u_n(:,end));
                    
                end
                
                
                for n = 1:length(obj.R)

                    [s, time] = SDESolver(problem, data, obj.N(end), 0, 1, obj.method{m}, obj.R(n), obj.M);
                    
                    obj.t              = [obj.t time];
                    
                    if obj.NoExact
                        
                        obj.Error.u_n(m, n) = mean(s.Solution.u_n(:,end));
                    
                        if n > 1

                            obj.Error.strong(m, n - 1) = abs(obj.Error.u_n(m,1) - obj.Error.u_n(m,n));
                        end
                        
                        if n > 2
                            
                            obj.Order.strong(m, n - 2) = log2(obj.Error.strong(m, n-2)/obj.Error.strong(m,n-1))/log2(obj.R(n - 1)/obj.R(n));
                            
                        end
                        
                        
                    else
                        
                        obj.Error.strong(m,n) = s.Solution.error_strong(:,end);
                        
                        if n > 1

                            obj.Order.strong(m, n - 1) = log2(obj.Error.strong(m, n-1)/obj.Error.strong(m,n))/log2(obj.R(n-1)/obj.R(n));

                        end
                    end
                    
                    clc
                    Table(obj, 'Error strong')
                    Table(obj, 'Order strong')
                    
                end
                
                dt_s = obj.R * 1./obj.N;
            
                subplot(length(obj.method), 1, m)
                loglog(dt_s, obj.Error.strong(m, :), dt_s, dt_s.^obj.Acc(m))
                legend(obj.method{m}, num2str(obj.Acc(m)))
                title(['Strong Convergence: ', obj.method{m}])

            end
                        
        end
        function obj = WeakConvergence(obj, problem, data)
            
            figure
            
            for m = 1:length(obj.method)
                
                if obj.NoExact
                    
                    [s, ~] = SDESolver(problem, data, obj.N(end), 0, 1, obj.method{m}, obj.R(1), obj.M);
                    
                    obj.Error.u_n(m, length(obj.N)) = mean(s.Solution.u_n(:, end));
                    
                end 
                
                for n = 1:length(obj.N) - obj.NoExact

                    [s, time] = SDESolver(problem, data, obj.N(n), 0, 1, obj.method{m}, obj.R(1), obj.M);

                    obj.t = [obj.t time];
                    
                    if obj.NoExact
                        
                        obj.Error.weak(m, n) = abs(mean(s.Solution.u_n(:, end)) - obj.Error.u_n(m, length(obj.N)));
                        
                        if n > 1
                            
                            obj.Order.weak(m, n - 1) = log2(obj.Error.weak(m, n-1)/obj.Error.weak(m,n)) / log2(obj.N(n)/obj.N(n - 1));
                            
                        end
                        
                    else
                        
                        obj.Error.weak(m,n) = s.Solution.error_weak(:,end);
                        
                        if n > 1

                            obj.Order.weak(m, n - 1) = log2(obj.Error.weak(m, n-1)/obj.Error.weak(m,n))/log2(obj.N(n)/obj.N(n - 1));

                        end
                    end
                    
                    clc
                    Table(obj, 'Error weak')
                    Table(obj, 'Order weak')
                    
                end
                
                dt_w = obj.R * 1./obj.N(1:end-obj.NoExact);
            
                
                subplot(length(obj.method), 1, m)
                loglog(dt_w, obj.Error.weak(m, :), dt_w, dt_w.^obj.Acc(m))
                legend(obj.method{m}, num2str(obj.Acc(m)))
                title(['Weak Convergence: ', obj.method{m}])
            end
            
            
            
        end
    end
end