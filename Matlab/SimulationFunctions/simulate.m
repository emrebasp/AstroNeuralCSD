
function result = simulate(modelp, simp, tSpan, initConditions, modeSpeed)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % "SIMULATE" FUNCTION
    %
    % RESULT = SIMULATE(PARAMETERS) provides the
    % result of the simulations of the astro-neural model. The result
    % correspond to the model and simulations parameters provided by MODELP
    % and SIMP, respectively. The result is a structur with
    %
    % RESULT.X : spatial grid
    % RESULT.T: time samples to be plotted
    % RESULT.Z: matrix containing the time integration result of the state
    % variables va, vi, ra, cK.
    %
    % This code was designed and edited by
    % Daniele Avitabile, VU Amsterdan, Netherlands
    % d.avitabile@vu.nl
    % Emre Baspinar, Inria, MathNeuro Team, France
    % emre.baspinar@inria.fr 
    % 
    % It was contributed by F. Campillo.
    %
    % Contact: emre.baspinar@inria.fr
    %    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Dependencies
    % 
    %
    % Simulation functions:
    %
    % Function LinearOperators
    % Function NeuralFieldWithDiffusible
    % Function TimeOutput
    % Function PopulationTransferFcn
    %
    %
    % Plot functions:
    %
    % Function PlotSolution
    % Function PlotHistory
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Initialization
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % clear all, close all, clc;
    % warning('off', 'MATLAB:dispatcher:UnresolvedFunctionHandle') % Turn off this warning,
    % % it does not have any influence on the code functionality.
    % 
    % % Add dependencies
    % addpath('Parameters/modelParameters');
    % addpath('Parameters/simulationParameters');
    % addpath('PlotFunctions');
    % addpath('SimulationFunctions');
    
    %==========================================================================
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Parameters
    % 
    % Model parameters (modelp_control) are for the model equations. 
    % Simulation parameters (simp) are for the discretization and the time integration
    % of the model. 
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Assign the simulation parameters
    nx = simp.nx;
    Lx = simp.Lx;
    % T = simp.T;
    % tStep = simp.tStep;
    % tSpan = [0:tStep:T]; % Output is sampled at the instants in tSpan.
    
    waveFrontThreshold = simp.waveFrontThreshold;
    errorVal = simp.errorVal;
    
    %==========================================================================

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Vector construction for the state variables
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Construct the vectors containing the discretized ve, vi, ra and cK. 
    iVe = [1:nx]'; % Indexes for ve vector (exc. population)
    iVi = nx+iVe;  % Indexes for vi vector (inh. population)
    iRa = nx+iVi;  % Indexes for ra vector (astro. population)
    iCK = nx+iRa;  % Indexes for cK vector (extracellular potassium concentration)
    
    idx = [iVe iVi iRa iCK]; % We put all indexes in idx, which has size(idx)=(nx, 4).    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Discretization in space, connectivity kernel and the finite differences 
    % 
    % x: position of nodes of the discretization of the cortical layer
    % kappa: connectivity kernel
    % Dxx: Second-order central finite difference to approximate the Laplacian
    % in the potassium equation.
    % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    convDomain = 'fourier'; % uncomment for convolutions in Fourier domain
    % convDomain = 'spatial'; % uncomment for convolutions in spatial domain
    [x, kappa, Dxx] = LinearOperators(nx, Lx, convDomain);
  
    %==========================================================================
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 
    % Simulation 
    %  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    
    % Initial conditions of the state variables
    z0 = initConditions;
    
    % Solve the system and plot its time evolution
    
    % Create a parent for the figures of time evolutions of ve, w and se
    
    % Define handle to right-hand side and time output function
    prob     = @(t,z) NeuralFieldWithDiffusible(t, z, modelp, Dxx, kappa, x, Lx, idx, convDomain);  % Time evolution of ve, vi, ra, cK  at the sampled instants by tSpan.
    timeOut = @(t,z,flag) TimeOutput(t, z, flag, true, x, modelp, simp, parentSolution, idx);       % Time evolution of ve, cK and Se for plots
    
    
    % Time integration of the model via ode45
    opts = odeset('OutputFcn', []);   
    [t,Z] = ode45(prob, tSpan, z0, opts);     % t gives the sampled instants in tSpan, 
                                              % Z is a (4 X nx) vector which contains
                                              % ve, vi, ra and cK at the sampled
                                              % instants.
    
    %--------------------------------------------------------------------------

    % Space-time plots
    
    % Assign the parameters for the transfer functions Se, Si, Sa
    theta1_e = modelp.theta1_e;
    theta2_e = modelp.theta2_e;
    theta3_e = modelp.theta3_e;
    theta1_i = modelp.theta1_i;
    theta2_i = modelp.theta2_i;
    theta3_i = modelp.theta3_i;
    theta1_a = modelp.theta1_a;
    theta2_a = modelp.theta2_a;
    theta3_a = modelp.theta3_a;
    theta4   = modelp.theta4;
    
    % Plot
    Se  = PopulationTransferFcn(Z(:,iVe), Z(:,iCK), theta1_e, theta2_e, theta3_e, theta4);   % Compute Se
    Si  = PopulationTransferFcn(Z(:,iVi), Z(:,iCK), theta1_i, theta2_i, theta3_i, theta4);   % Compute Si
    Saa = PopulationTransferFcn(Z(:,iRa), Z(:,iCK), theta1_a, theta2_a, theta3_a, theta4);   % Compute Saa
    
    %--------------------------------------------------------------------------
    
    if modeSpeed == 1
        % Compute the propagation speed 
        
        % We use the exc. population activity to compute the propagation speed.
        Cx = simp.Cx; % Scale factor to fit the model propagation speed to the experimentally obtained speed
        hx = Lx/(nx-1); % Distance between the node of the discretization of cortical layer
        dist = abs(Z(end,idx(:, 1)) - waveFrontThreshold); % Difference between ve and the threshold value 
                                                           % at every node, and at the end of the simulation 
        minDist = min(dist); % Minimum of the differences.
        tKCl = modelp.tKCl; % Instant when the potassium puff is started to be applied.
        
        if minDist>errorVal  % Check if this min is below the tolerance.
            disp('No propagating excitatory wave! You might try to increase errorVal.')
        else
            indFind = find(dist == minDist) % The node index corresponding to the min of the differences.
            % [~, waveFrontNode] = min(dist);
            waveFrontNode = max(indFind); % Propagation is symmetric around the center. Pick the node propagating to the right on the cortical layer.
            waveFrontDistance = abs(-Lx/2 + (waveFrontNode-1)*hx); % Distance between the center and detected wavefront
            propagationSpeed = Cx * (waveFrontDistance/(tSpan(end)-tKCl)); % Propagation speed
            fprintf('Propagation speed: %.2f mm/min\n\n', propagationSpeed); % Display the propagation speed
        end

        result.propagationSpeed = propagationSpeed;

    end
    %--------------------------------------------------------------------------

    result.x = x;
    result.t = t;
    result.Z = Z;
    result.Se = Se;
    result.Si = Si;
    result.Saa = Saa;

end


