classdef floris<handle
    properties
        inputData
        outputData
        outputDataAEP
    end
    methods        
        %% Initialization
        function self = init(self,modelType,turbType,siteType)
            addpath('functions');
            
            % Default setup
            if exist('siteType') == 0;  siteType  = '9turb';   end;
            if exist('turbType') == 0;  turbType  = 'NREL5MW'; end;
            if exist('modelType')== 0;  modelType = 'default'; end;
            
            % Call function
            self.inputData = floris_loadSettings(modelType,turbType,siteType);
        end
       
       
        
        %% FLORIS single execution
        function [self,outputData] = run(self,inputData)
            % Check if inputData specified manually. If not, use internal inputData.
            if exist('inputData') == 0; inputData = self.inputData; end;
            
            % Check if init() has been run at least once before.
            if isstruct(inputData)== 0;
                disp([' Please initialize the model before simulation by the init() command.']);
                return; end;
            
            % Run FLORIS simulation
            [self.inputData,self.outputData] = floris_run(inputData);
            
            % Results saved internally, but also returns externally if desired.
            if nargout > 0; outputData = self.outputData; end;
        end
        
        
        %% Optimization functions
        function [axialIndOpt,self] = optimizeAxInd(self,inputData);
            % Check if inputData specified manually. If not, use internal inputData.
            if exist('inputData') == 0; inputData = self.inputData; end;            
            if isstruct(inputData)== 0; % Check if init() has been run at least once before.
                disp([' Please initialize the model before simulation by the init() command.']);
                return; end;
            
            % Cost function
            function J = costFunction(axialIndIn,inputData)
                inputData.axialInd  = axialIndIn;
                [~,outputData]      = floris_run(inputData,0);
                J                   = -sum(outputData.power);
            end
            cost  = @(x)costFunction(x,self.inputData);
            lb    = 0.000*ones(length(self.inputData.axialInd),1); 
            ub    = 0.333*ones(length(self.inputData.axialInd),1);
            
            % Optimizer settings and optimization execution
            %options = optimset('Display','final','MaxFunEvals',1000 ); % Display nothing
            options = optimset('Algorithm','sqp','Display','final','MaxFunEvals',1000,'PlotFcns',{@optimplotx, @optimplotfval} ); % Display convergence
            axialIndOpt = fmincon(cost,self.inputData.axialInd,[],[],[],[],lb,ub,[],options);
                      
            % Simulated annealing
            %options = optimset('Display','iter','MaxFunEvals',1000,'PlotFcns',{@optimplotx, @optimplotfval} ); % Display convergence
            %axialIndOpt = simulannealbnd(cost,self.inputData.axialInd,lb,ub,options);
            
            % Overwrite current settings with optimized
            self.inputData.axialInd = axialIndOpt;
        end;
        
        
        function [yawAnglesOpt,self] = optimizeYaw(self,inputData);
            % Check if inputData specified manually. If not, use internal inputData.
            if exist('inputData') == 0; inputData = self.inputData; end;            
            if isstruct(inputData)== 0; % Check if init() has been run at least once before.
                disp([' Please initialize the model before simulation by the init() command.']);
                return; end;
            
            % Cost function
            function J = costFunction(yawAnglesIn,inputData)
                inputData.yawAngles = yawAnglesIn;
                [~,outputData]      = floris_run(inputData,0);
                J                   = -sum(outputData.power);
            end
            cost         = @(x)costFunction(x,self.inputData);
            
            % Optimizer settings and optimization execution
            %options = optimset('Display','final','MaxFunEvals',1000 ); % Display nothing
            options = optimset('Display','final','MaxFunEvals',1000,'PlotFcns',{@optimplotx, @optimplotfval} ); % Display convergence
            
            % Yaw angles are constrained between -25 and +25 degrees
            lb = deg2rad(-25)*ones(length(self.inputData.yawAngles),1); 
            ub = deg2rad(+25)*ones(length(self.inputData.yawAngles),1);
            yawAnglesOpt = fmincon(cost,self.inputData.yawAngles,[],[],[],[],lb,ub,[],options);
            
            % Overwrite current settings with optimized
            self.inputData.yawAngles = yawAnglesOpt;
        end;
        
        
        %% Visualize single FLORIS simulation results
        function [] = visualize(self,plotLayout,plot2D,plot3D)
            inputData  = self.inputData;
            outputData = self.outputData;
            
            if isstruct(outputData) == 0
                disp([' outputData is not (yet) available/not formatted properly.' ...
                      ' Please run a (single) simulation, then call this function.']);
                return;
            end; 
            
            % Default visualization settings, if not specified
            if exist('plotLayout') == 0; plotLayout = true;  end;
            if exist('plot2D')     == 0; plot2D     = true;  end;
            if exist('plot3D')     == 0; plot3D     = false; end;
            
            % Set visualization settings
            inputData.plotLayout      = plotLayout;
            inputData.plot2DFlowfield = plot2D;
            inputData.plot3DFlowfield = plot3D;
            
            floris_visualization(inputData,outputData);
        end;
        
        
        %% Run FLORIS AEP calculations (multiple wind speeds and directions)
        function [self,outputDataAEP] = AEP(self,windRose)
            % WindRose is an N x 2 matrix with uIf in 1st column and 
            % vIf in 2nd. The simulation will simulate FLORIS for each row.
            inputData = self.inputData;
            
            % Simulate over each uIf-vIf set (matrix row)
            for i = 1:size(windRose,1);
                inputData.uInfIf      = WS_range(windRose,1);
                inputData.vInfIf      = WS_range(windRose,2);
                self.outputDataAEP{i} = floris_run(inputData);
            end;
            
            % Results saved internally, but also returns externally if desired.
            if nargout > 0; outputDataAEP = self.outputDataAEP; end;
        end        
    end
end