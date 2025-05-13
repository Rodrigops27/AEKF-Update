%% Implementing AEKF Lesson 3.3.6
% The original source code was developed by Prof. Plett, UCCS
% for the Applied Kalman Filtering Specialization.
% It is reproduced here with proper attribution for academic and
% educational use.

% --------------------------------------------------------------------
% Load data from simulation of system dynamic model.
% Data file contains the variables: model dT t u x z
% Note that this is a longer simulation of the model than we have used before.
% The longer simulation gives time to adapt SigmaV and SigmaW
clear; % 1-line Update to ensure clean start.
load readonly/simOutLong.mat
rms = @(x) sqrt(mean((x).^2)); % define root-mean-squared function

% Determine number of states and timesteps
[nx,nt]    = size(x);
% Reserve storage for results
xhatstore  = zeros(nx,nt);
boundstore = zeros(nx,nt);
Qstore = zeros(nx^2,nt);
Rstore = zeros(1,nt);

% Covariance values
SigmaX0 = 0.1*eye(nx,nx); % uncertainty of initial state
xhat0   = zeros(nx,1);

% Create aekfData structure and initialize variables
aekfData = initAEKF(xhat0,SigmaX0,u(1),model);

% Now, enter loop for remainder of time, where we update the EKF
% once per sample interval
for k = 2:length(t)
    uk = u(k); % "measure" input force
    zk = z(k); % "measure" nonlinear output

    % Compute estimate of present state vector and its bounds
    [xhatstore(:,k),boundstore(:,k),Qstore(:,k),Rstore(k),aekfData] = ...
        iterAEKF(uk,zk,dT,aekfData);
    % if mod(aekfData.iter,10000) == 0
    %   fprintf('Simulated %d out of %d iterations\n',aekfData.iter,nt);
    % end
end

% Isolate the final 1000 samples for plotting
tseg = t(end-1000:end); tseg = tseg - tseg(1);
xseg = x(:,end-1000:end);
xhatseg = xhatstore(:,end-1000:end);
boundseg = boundstore(:,end-1000:end);
errseg = xseg - xhatseg;

% Compute and output some statistics
fprintf('RMS state estimation error = %g\n',rms(errseg(:)));
ind = find(abs(errseg(:))>boundseg(:));
fprintf('Time error outside bounds = %g%%\n',length(ind)/length(boundseg(:))*100);

%% Plot the states with confidence bounds
% 1-line Update:
plotting; % refactored for readiblity

%% AEKF Functions
% AEKF initialization function, called once every simulation
function aekfData = initAEKF(xhat0,SigmaX0,priorU,model)
aekfData.xhat   = xhat0;        % Initial state description
aekfData.SigmaX = SigmaX0;      % Initial estimation-error covariance
aekfData.SigmaV = model.SigmaV; % Sensor-noise covariance
aekfData.SigmaW = model.SigmaW; % Process-noise covariance
aekfData.priorU = priorU;       % Previous value of input force
aekfData.model  = model; % Store model data structure too

%% For adaptive code: Storage of innovations and residues
NV = 1000; % Length of buffer of past SigmaV to average
aekfData.Vstore = zeros(1,NV);
NW = 1000; % Length of buffer of past SigmaW to average
aekfData.Wstore = zeros(length(xhat0)^2,NW);

aekfData.iter = 0; % The present iteration number "k"
end

% AEKF iteration function, called once every measurement interval
function [xhat,bounds,SigmaW,SigmaV,aekfData] = iterAEKF(uk,zk,dT,aekfData)
model = aekfData.model;

% Get data stored in ekfData structure
priorU = aekfData.priorU;
SigmaX = aekfData.SigmaX;
SigmaV = aekfData.SigmaV;
SigmaW = aekfData.SigmaW;
xhat   = aekfData.xhat;

% EKF Step 0: Compute Ahat[k-1]
Ahat = [1 dT; -model.k1*dT/model.m-3*model.k2*dT*xhat(1)^2/model.m ...
    1-model.b*dT/model.m];

% Step 1a: State estimate time update
xhat = [1 dT; -model.k1*dT/model.m 1-model.b*dT/model.m]*xhat + ...
    [0; -model.k2*dT*xhat(1)^3/model.m] + [0; dT/model.m]*priorU;

% 1-line Update:
Sx_prev = SigmaX; % To use prev step.
% Step 1b: Error covariance time update
%          sigmaminus(k) = Ahat(k-1)*sigmaplus(k-1)*Ahat(k-1)' + ...
%                          Bhat(k-1)*sigmawtilde*Bhat(k-1)'
% Note that Bhat is taken into account already in SigmaW by initEKF.m
SigmaX = Ahat*SigmaX*Ahat' + SigmaW;

% Step 1c: Output estimate
zhat = atan(xhat(1)/model.d);

% Step 2a: Estimator gain matrix
Chat = [model.d/(model.d^2 + xhat(1)^2) 0]; Dhat = 1;
SigmaZ = Chat*SigmaX*Chat' + Dhat*SigmaV*Dhat';
L = SigmaX*Chat'/SigmaZ;

% Step 2b: State estimate measurement update
mu = zk - zhat; % innovation.  Use to check for sensor errors...
xhat = xhat + L*mu;
zhatplus = atan(xhat(1)/model.d);
r = zk - zhatplus; % residue

% Step 2c: Error covariance measurement update
SigmaX = SigmaX - L*SigmaZ*L';
[~,S,V] = svd(SigmaX); % Implement Higham method to help robustness
HH = V*S*V';
SigmaX = (SigmaX + SigmaX' + HH + HH')/4;

% Update:
%% Testing adaptation techniques
alpha = 0.99; % simple filter pole location

aekfData.mode = 3;
% 1: Original (sum), 2: Full def. of SigmaW, 3: Full def. +Higham
switch aekfData.mode
    case 1 % Original
        % SigmaW: cumsum + filter
        newW = (L*mu)*(L*mu)';
        aekfData.Wstore = [aekfData.Wstore(:,2:end), newW(:)]; % Buffer result
        if aekfData.iter > size(aekfData.Wstore,2) % Smooth with buffer AND filter
            newW = sum(aekfData.Wstore,2); % WARNING: Should be mean and not sum
            aekfData.SigmaW(:) = alpha*aekfData.SigmaW(:) + (1-alpha)*newW;
        end
        % For adapting SigmaV
        newV = r*r' + Chat*SigmaX*Chat';
        aekfData.Vstore = [aekfData.Vstore(2:end), newV(:)]; % Buffer result
        if aekfData.iter > size(aekfData.Vstore,2) % Smooth with buffer AND filter
            newV = mean(aekfData.Vstore,2);
            aekfData.SigmaV(:) = alpha*aekfData.SigmaV(:) + (1-alpha)*newV;
        end

    case 2 % GenFilt + Std SigmaV
        % SigmaW: MAvg + filter
        newW = (L*mu)*(L*mu)' + (SigmaX - Ahat*Sx_prev*Ahat'); % Generalized
        aekfData.Wstore = [aekfData.Wstore(:, 2:end), newW(:)];
        if aekfData.iter > size(aekfData.Wstore, 2)
            newW = mean(aekfData.Wstore, 2);
            aekfData.SigmaW(:) = alpha * aekfData.SigmaW(:) + (1 - alpha) * newW;
        end
        % Std adapting SigmaV
        newV = r*r' + Chat*SigmaX*Chat';
        aekfData.Vstore = [aekfData.Vstore(2:end), newV(:)]; % Buffer result
        if aekfData.iter > size(aekfData.Vstore,2) % Smooth with buffer AND filter
            newV = mean(aekfData.Vstore,2);
            aekfData.SigmaV(:) = alpha*aekfData.SigmaV(:) + (1-alpha)*newV;
        end

    case 3 % GenFilt + Std SigmaV
        % SigmaW: Higham + MAvg + filter
        newW = (L*mu)*(L*mu)' + (SigmaX - Ahat*Sx_prev*Ahat'); % Generalized
        % Robustifying
        [~,S,V] = svd(newW); % Implement Higham method to help robustness
        HH = V*S*V';
        newW = (newW + newW' + HH + HH')/4;

        aekfData.Wstore = [aekfData.Wstore(:, 2:end), newW(:)];
        if aekfData.iter > size(aekfData.Wstore, 2)
            newW = mean(aekfData.Wstore, 2);
            aekfData.SigmaW(:) = alpha * aekfData.SigmaW(:) + (1 - alpha) * newW;
        end
        % Std adapting SigmaV
        newV = r*r' + Chat*SigmaX*Chat';
        aekfData.Vstore = [aekfData.Vstore(2:end), newV(:)]; % Buffer result
        if aekfData.iter > size(aekfData.Vstore,2) % Smooth with buffer AND filter
            newV = mean(aekfData.Vstore,2);
            aekfData.SigmaV(:) = alpha*aekfData.SigmaV(:) + (1-alpha)*newV;
        end
end
% end of Update.

% Save data in ekfData structure for next time...
aekfData.priorU = uk;
aekfData.SigmaX = SigmaX;
aekfData.xhat = xhat;
bounds = 3*sqrt(diag(SigmaX));

aekfData.iter = aekfData.iter + 1;
SigmaW = aekfData.SigmaW(:);
SigmaV = aekfData.SigmaV;
end
