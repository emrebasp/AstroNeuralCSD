function Dt = NeuralFieldWithDiffusible(t, z, modelp, Dxx, kappa, x, Lx, idx, check)

% Dt = NEURALFIELDWITHDIFFUSIBLE(T, Z, MODELP, DXX, KAPPA, X, LX, IDX) performs
% the time integration of the model via ode45. It provides the numerical
% result F of the state variables ve, vi, ra, cK at each instant given by T.
% F is a (4 X nx) dimensional vector, where nx is the number of discretization
% nodes on the cortical layer. F has the following form:
%     
%         |  dve(x_1,t)/dt   |
%         |      .           |
%         |      .           |
%         |  dve(x_nx,t)/dt  |
%         |  dvi(x_1,t)/dt   |   
%         |      .           |
%         |      .           |
%         |  dvi(x_nx,t)/dt  |
% Dt(t) = |  dra(x_1,t)/dt   | .
%         |      .           |
%         |      .           |
%         |  dra(x_nx,t)/dt  |
%         |  dcK(x_1,t)/dt   |
%         |      .           |
%         |      .           |
%         |  dcK(x_nx,t)/dt  |
%          
% 
% NEURALFIELDWITHDIFFUSIBLE(T, Z, MODELP, DXX, KAPPA, X, LX, IDX): T is the time span
% containing the instants at which the state variables ve, vi, va, cK will
% be registered. Z is a (4 X nx) dimensional vector which contains all state variables
% discretized on the cortical layer. It has the following form:
%         
%        |  ve(x_1,t)  |
%        |      .      |
%        |      .      |
%        |  ve(x_nx,t) |
%        |  vi(x_1,t)  |   
%        |      .      |
%        |      .      |
%        |  vi(x_nx,t) |
% Z(t) = |  ra(x_1,t)  | .
%        |      .      |
%        |      .      |
%        |  ra(x_nx,t) |
%        |  cK(x_1,t)  |
%        |      .      |
%        |      .      |
%        |  cK(x_nx,t) |   
% 
% MODELP is a structure which contains all the model parameters. DXX is a 
% (nx X nx) dimensional matrix which represents the
% Laplacian approximated via second-order central finite differences. It
% has the following form:
%                                      
%       | -2   1   0   0   0  . . .   0   0   1 |
%       |  1  -2   1   0   0  . . .   0   0   0 |
%       |  0   1  -2   1   0  . . .   0   0   0 |
%       |  0   0   0   .                  .   . |
%       |  .   .   .       .              .   . |
% Dxx = |  .   .   .           .          .   . | .
%       |  .   .   .                .     .   . |
%       |  .   .   .                    .       |
%       |  0   0   0   0   0  . . .   1  -2   1 |
%       |  1   0   0   0   0  . . .   0   1  -2 |
%  
% KAPPA is the Fourier transform of the connectivity kernel. X is
% a nx dimensional vector providing the discretization nodes on the
% cortical layer. Lx is the length of the cortical layer.
% IDX is a (4 X nx) vector providing the indexes of the state variables
% ve, vi, ra, cK corresponding to the discretizatoin nodes in X.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Dependencies
%
% PopulationTransferFc
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Assign the parameters
nx = size(z,1)/4;

% Split state variables
iVe = idx(:,1);  % indexes of discretized ve (exc. population)
ve  = z(iVe);    % discretized ve
iVi = idx(:,2);  % indexes of discretized vi (inh. population)
vi  = z(iVi);    % discretized vi
iRa = idx(:,3);  % indexes of discretized ra (astro. population)
ra  = z(iRa);    % discretized ra
iCK = idx(:,4);  % indexes of discretized cK (extracellular potassium concentration)  
cK  = z(iCK);     % discretized cK

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Convolutions between connectivity kernel \kappa and 
% population transfer functions Se, Si and Sa
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of the population transfer functions
theta4          = modelp.theta4;       % threshold of the transfer functions Se, Si, Sa
theta3_e        = modelp.theta3_e;     % sharpness of nonlinearity of Se 
theta3_i        = modelp.theta3_i;     % sharpness of nonlinearity of Si
theta3_a        = modelp.theta3_a;     % sharpness of nonlinearity of Sa
theta1_e        = modelp.theta1_e;     % initial potassium level of hyperactivation of the exc. population
theta2_e        = modelp.theta2_e;     % initial potassium level of depolarization block of the exc. population
theta1_i        = modelp.theta1_i;     % initial potassium level of hyperactivation of the inh. population
theta2_i        = modelp.theta2_i;     % initial potassium level of depolarization block of the inh. population
theta1_a        = modelp.theta1_a;     % initial potassium level of increased potassium clearance of the astro. population
theta2_a        = modelp.theta2_a;     % initial potassium level of no potassium clearance

%--------------------------------------------------------------------------

if check == 'fourier'
    
    % Convolutions are via Fourier transform. Here ifftshift is needed 
    % due to Matlab convention.
    Se      = PopulationTransferFcn(ve, cK, theta1_e, theta2_e, theta3_e, theta4);     % Output of exc. population transfer fcn. Se
    SeConv  = Lx/nx * ifftshift(real(ifft(fft(Se) .* kappa)));                         % Convolution with connectivity kernel \omega, wHat is fft(\omega)
    Si      = PopulationTransferFcn(vi, cK, theta1_i, theta2_i, theta3_i, theta4);     % Output of inh. population transfer fcn. Si
    SiConv  = Lx/nx * ifftshift(real(ifft(fft(Si) .* kappa)));                         % Convolution with the connectivity kernel
    San     = PopulationTransferFcn(ve+vi, cK, theta1_a, theta2_a, theta3_a, theta4);  % Output of astrocyte population transfer fcn. Sa in San (astro-neuron interactions)
    SanConv = Lx/nx * ifftshift(real(ifft(fft(San) .* kappa)));                        % Convolution with the connectivity kernel
    Saa     = PopulationTransferFcn(ra, cK, theta1_a, theta2_a, theta3_a, theta4);     % Output of astrocyte population transfer fcn. Sa in Saa (astro-astro interactions)
    SaaConv = Lx/nx * ifftshift(real(ifft(fft(Saa) .* kappa)));                        % Convolution with the connectivity kernel

elseif check == 'spatial'
    
    % Convolutions in spatial domain
    Se      = PopulationTransferFcn(ve, cK, theta1_e, theta2_e, theta3_e, theta4);     % Output of exc. population transfer fcn. Se
    SeConv  = Lx/nx * conv(Se, wConnectivity,'same');                                  % Convolution with connectivity kernel \omega, wHat is fft(\omega)
    Si      = PopulationTransferFcn(vi, cK, theta1_i, theta2_i, theta3_i, theta4);     % Output of inh. population transfer fcn. Si
    SiConv  = Lx/nx * conv(Si, wConnectivity,'same');                                  % Convolution with the connectivity kernel
    San     = PopulationTransferFcn(ve+vi, cK, theta1_a, theta2_a, theta3_a, theta4);  % Output of astrocyte population transfer fcn. Sa in San (astro-neuron interactions)
    SanConv = Lx/nx * conv(San, wConnectivity,'same');                                 % Convolution with the connectivity kernel
    Saa     = PopulationTransferFcn(ra, cK, theta1_a, theta2_a, theta3_a, theta4);     % Output of astrocyte population transfer fcn. Sa in Saa (astro-astro interactions)
    SaaConv = Lx/nx * conv(Saa, wConnectivity,'same');                                 % Convolution with the connectivity kernel

end

%========================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potassium-dependent activation function F: extracellular potassium effects on cell
% population activity represented by ve, vi, ra
%
% Inputs F
% cK: extracellular potassium concentration
% alpha1, alpha2: parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of F
alpha1_e      = modelp.alpha1_e;     % sharpness of the nonlinerarity of Fe in ve equation
alpha1_i      = modelp.alpha1_i;     % sharpness of the nonlinerarity of Fi in vi equation
alpha1_a      = modelp.alpha1_a;     % sharpness of the nonlinerarity of Fa in va equation
alpha2_e      = modelp.alpha2_e;     % potassium threshold of Fe appearing in ve equation
alpha2_i      = modelp.alpha2_i;     % potassium threshold of Fi appearing in vi equation
alpha2_a      = modelp.alpha2_a;     % potassium threshold of Fa function appearing in ra equation

% Weights of F functions
gamma_eK = modelp.gamma_eK; % weight for Fe
gamma_iK = modelp.gamma_iK; % weight for Fi
gamma_aK = modelp.gamma_aK; % weight for Fa

% Definition of F
F = @(cK, alpha1, alpha2) 1 ./ (1 + exp(- alpha1 * (cK - alpha2) )); 

% Define Fe, Fi, Fa
Fe = @(cK) F(cK, alpha1_e, alpha2_e);
Fi = @(cK) F(cK, alpha1_i, alpha2_i);
Fa = @(cK) F(cK, alpha1_a, alpha2_a);

% Weights of potassium-dependent activation functions Fe, Fi, Fa
gamma_eK = modelp.gamma_eK;         % weight for Fe
gamma_iK = modelp.gamma_iK;         % weight for Fi
gamma_aK = modelp.gamma_aK;         % weight for Fa

%======================================================================== 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Potassium source function rho: contribution of neural spiking and astrocyte activity
% to the increase in extracellular potassium concentration cK
%
% Inputs rho 
% S: output of the population transfer function
% beta1, beta2: parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of rho
beta1_e = modelp.beta1_e;     % sharpness of the nonlinearity of rho function weighted by gamma_e in potassium equation
beta2_e = modelp.beta2_e;     % transfer function threshold for rho function weighted by gamma_e in potassium equation
beta1_i = modelp.beta1_i;     % sharpness of the nonlinearity of rho function weighted by gamma_i in potassium equation
beta2_i = modelp.beta2_i;     % transfer function threshold for rho function weighted by gamma_i in potassium equation
beta1_a = modelp.beta1_a;     % sharpness of the nonlinearity of rho function weighted by gamma_a in potassium equation
beta2_a = modelp.beta2_a;     % transfer function threshold for rho function weighted by gamma_a in potassium equation  

% Definition of rho
rho = @(S, beta1, beta2) 1 ./ cosh(beta1*(S-beta2));

%========================================================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% KCl stimulus input I: potassium puff which triggers CSD
%
% Inputs I
% x: position on the cortical layer
% t: time
% eta1, eta2, tKCl: parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of the external stimulus I
eta1           = modelp.eta1;  % amplitude of the external stimulus I
eta2           = modelp.eta2;  % sharpness of the nonlinearity (in space) of the external stimulus
tKCl           = modelp.tKCl; % instant when the external stimulus I starts

% Definition of I
I = @(x, t, eta1, eta2, tKCl) exp(-0.5.*(t-tKCl)).*eta1*(t>=tKCl)./cosh(eta2 * x).^2;    

%========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Model equations: time step for the time integration of the model
%
% The system equations are discretized at nx points on the one-dimensional
% cortical layer. Therefore, it can be written as a 4nx dimensional system 
% of ODEs here.
% 
% Dt(iVe): nx-dimensional vector which contains the dve/dt at each
% discretization point on the one-dimensional cortical layer.
% 
% Dt(iVi): nx-dimensional vector which contains the dvi/dt at each
% discretization point on the one-dimensional cortical layer.
% 
% Dt(iRa): nx-dimensional vector which contains the dra/dt at each
% discretization point on the one-dimensional cortical layer.
% 
% Dt(iCK): nx-dimensional vector which contains the dcK/dt at each
% discretization point on the one-dimensional cortical layer.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters for the connectivity and extracellular potassium
% accumulation
gamma_e          = modelp.gamma_e;  % gamma_e: contribution weight of exc. population spiking to extracellular potassium accumulation
gamma_i          = modelp.gamma_i;  % gamma_i: contribution weight of inh. population spiking to extracellular potassium accumulation
gamma_a          = modelp.gamma_a;  % gamma_a: contribution weight of astrocyte potassium clearance from the extracellular matrix

gamma_ei         = modelp.gamma_ei; % gamma_ei: connection weight from inh. to exc. neurons
gamma_ie         = modelp.gamma_ie; % gamma_ie: connection weight from exc. to inh. neurons
gamma_ee         = modelp.gamma_ee; % gamma_ee: connection weight from exc. to exc. neurons
gamma_ii         = modelp.gamma_ii; % gamma_ii: connection weight from inh. to inh. neurons
gamma_an         = modelp.gamma_an; % gamma_an: connection weight from all neurons to astrocytes
gamma_aa         = modelp.gamma_aa; % gamma_aa: connection weight from astrocytes to astrocytes

% Initialization
Dt    = zeros(size(z));
sigma = modelp.sigma;   % diffusion constant

amp = 1;
cK_th = 1.5;
offset = amp * cK_th;

% Model equations
Dt(iVe) = -ve + gamma_ee*SeConv - gamma_ei*SiConv + gamma_eK*Fe(cK);
Dt(iVi) = -vi + gamma_ie*SeConv - gamma_ii*SiConv + gamma_iK*Fi(cK);
Dt(iRa) = -ra + gamma_an*SanConv + gamma_aa*SaaConv + gamma_aK*Fa(cK);
% Dt(iCK) = sigma*Dxx*cK + gamma_e*rho(Se, beta1_e, beta2_e) + gamma_i*rho(Si, beta1_i, beta2_i) - (2.5 - iCK.^1.1) .* (gamma_a*rho(Saa, beta1_a, beta2_a))...
%           + I(x, t, eta1, eta2, tKCl);
Dt(iCK) = sigma*Dxx*cK + gamma_e*rho(Se, beta1_e, beta2_e) + gamma_i*rho(Si, beta1_i, beta2_i) - (-amp * (z(iCK)-cK_th) .^2 + offset) .* (gamma_a*rho(Saa, beta1_a, beta2_a))...
          + I(x, t, eta1, eta2, tKCl);

end
 
