classdef error_obj
    % Error object
    properties
        
        N;     % Matrix of number of grid points
        Acc;   % Vector of order of accuracy
        t;     % Matrix of times
        Error; % Matrix of L_2 errors
        Order; % Matrix of computed order of accuracy
    end
    methods
        function obj = initialize(obj, problem, data)
            
            obj.t     = zeros(length(obj.Acc), size(obj.N, 2));
            obj.Error = zeros(length(obj.Acc), size(obj.N, 2));
            obj.Order = zeros(length(obj.Acc), size(obj.N, 2) - 1);
            format long
            
            for n = 1:size(obj.N, 2)
                for m = 1:length(obj.Acc)
                    tic;
                    %% Initializing grid
                    grid          = grid_obj;
                    grid.minv     = [0, 0, 0];
                    grid.maxv     = [1, 1, 1];
                    grid.N        = [obj.N(1, n) obj.N(2, n) 100];
                    grid.geometry = 'Circle Sector';
                    grid          = createGrid(grid);
                    
                    %% Initializing UQ
                    uq                         = uq_obj;
                    uq.Distribution.minv       = -1;
                    uq.Distribution.maxv       =  1;
                    uq.Grid.N                  = 0;
                    uq.Grid.dim                = 1;
                    uq.method                  = {'NI'};
                    uq.Distribution.name       = 'Uniform';
                    uq.Distribution.parameters = [0, 1];
                    uq                         = createGrid(uq);
                    
                    %% Initialize problem
                    
                    p = problem_obj;
                    p.type.problem = problem; 
                    p.type.data    = data;
                    p = initialize(p);
                    p = createMatrices(p, grid);
                    p = createData(p, grid);
                    
                    %% Initializing SBP operators
                    sbp     = sbp_obj;
                    sbp.acc = obj.Acc(m);
                    sbp     = createOperators(sbp, grid);
                    
                    %% Initializing scheme
                    s = scheme_obj;
                    s = createBoundaryOperator(s, grid, sbp, p);
                    s = createPenalties(s);
                    s = createAMatrix(s, grid, sbp, p, uq);
                    s = solveUsingRK4(s, grid, p, uq);
                    s = computeError(s, grid, sbp, p);
                    
                    %% Results: Statistics
                    uq.statistics.mean_a     = 0; %mean_a;
                    uq.statistics.variance_a = 0; %var_a;
                    [~, ~, ~] = functionals(uq, grid, s.Solution.u_ns);
                    obj.t(m,n) = toc;
                    obj.Error(m, n) = s.Solution.e_L2;
                    
                    if n > 1
                        
                        obj.Order(m, n-1) = log2(obj.Error(m, n-1)/obj.Error(m,n))/log2(obj.N(1, n)/obj.N(1, n-1));
                        
                    end
                    
                    clc
                    Table(obj, 'Error')
                    Table(obj, 'Order')
                    Table(obj, 'Time')
                    
                end
            end
        end
        function Table(obj, str)
            
            xVec = obj.N(1,:);
            order = obj.Acc;
            
            if strcmp(str, 'Order')
                
                titleH    = ' Order Of Accuracy ';
                spacing   = '\t';
                xVec(end) = [];
                data      = obj.Order;
                
            elseif strcmp(str, 'Error')
                
                titleH  = ' Error ';
                spacing = '\t\t';
                data    = obj.Error;
                
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
            fprintf('Order/(Nx,Ny)')
            fprintf('\t')
            
            for i = 1:length(xVec)
                
                fprintf([num2str(xVec(i)) ',' num2str(xVec(i))])
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
    end
end