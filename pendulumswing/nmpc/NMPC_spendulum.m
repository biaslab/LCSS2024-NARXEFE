%% Swing-up Control of a Pendulum Using Nonlinear Model Predictive Control
% Copyright 2016-2020 The MathWorks, Inc.

close all;
clear all;

%% Product Requirement
if ~mpcchecktoolboxinstalled('optim')
    disp('Optimization Toolbox is required to run this example.')
    return
end

%% Create Nonlinear MPC Controller
nx = 2;
ny = 1;
nu = 1;
nlobj = nlmpc(nx, ny, nu);

%% Sample time
Ts = 0.1;
nlobj.Ts = Ts;

%% Horizons
nlobj.PredictionHorizon = 5;
nlobj.ControlHorizon = 5;

%% Specify Nonlinear Plant Model
nlobj.Model.StateFcn = "pendulumDT0";
nlobj.Model.IsContinuousTime = false;
nlobj.Model.NumberOfParameters = 1;

%% Outputs
nlobj.Model.OutputFcn = 'pendulumOutputFcn';
nlobj.Jacobian.OutputFcn = @(x,u,Ts) [1 0];

%% Weights
nlobj.Weights.OutputVariables = [1];
nlobj.Weights.ManipulatedVariables = 1e-3;
nlobj.Weights.ManipulatedVariablesRate = 0.2;

%% Control limits
nlobj.MV.Min = -10;
nlobj.MV.Max = 10;

%% Validate Nonlinear MPC Controller
theta0 = 0.0;
theta_dot0 = 0.0;
x0 = [theta0; theta_dot0];
u0 = 0.0;
validateFcns(nlobj,x0,u0,[],{Ts});

%% State Estimation
EKF = extendedKalmanFilter(@pendulumStateFcn, @pendulumMeasurementFcn);

%% Closed-Loop Simulation in MATLAB(R)
x = [0.0; 0.0];
y = [x(1)];
EKF.State = x;
mv = 0;

%% Goal observation
yref = [pi];

%% Optimal control parameters
nloptions = nlmpcmoveopt;
nloptions.Parameters = {Ts};

%% Run the simulation for |20| seconds.
Duration = 20;
hbar = waitbar(0,'Simulation Progress');
x_NMPC = x;
y_NMPC = y;
u_NMPC = mv;
for ct = 1:(20/Ts)
    
    % Correct previous prediction using current measurement.
    xk = correct(EKF, y);

    % Compute optimal control moves.
    [mv,nloptions,info] = nlmpcmove(nlobj,xk,mv,yref,[],nloptions);

    % Predict prediction model states for the next iteration.
    predict(EKF, [mv; Ts]);

    % Implement first optimal control move and update plant states.
    x = pendulumDT0(x,mv,Ts);

    % Generate sensor data with some white noise.
    y = x(1) + randn(1)*1e-3;

    % Save plant states for display.
    y_NMPC = [y_NMPC y];
    x_NMPC = [x_NMPC x];
    u_NMPC = [u_NMPC mv];

    waitbar(ct*Ts/20,hbar);
end
close(hbar)

%% Plot the closed-loop response.
% figure
% subplot(1,2,1)
% plot(0:Ts:Duration,x_NMPC(1,:))
% xlabel('time')
% ylabel('theta')
% title('pendulum angle')
% subplot(1,2,2)
% plot(0:Ts:Duration,x_NMPC(2,:))
% xlabel('time')
% ylabel('thetadot')
% title('pendulum velocity')

figure
subplot(1,2,1)
plot(0:Ts:Duration,y_NMPC)
xlabel('time')
ylabel('theta')
subplot(1,2,2)
plot(0:Ts:Duration,u_NMPC)
xlabel('time')
ylabel('torque')

%% Store results
save('../results/NMPC.mat', 'x_NMPC', 'y_NMPC', 'u_NMPC'); 

