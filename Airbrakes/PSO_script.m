clc;
% clear;
close all;

% WARNING
% to avoid running into problems while running this script, a full run of
% the simulator before starting the tuning process is suggested

% INFO
% this PSO script is used for tuning the P, I, D and N values used in
% the PID which controls the air brakes control system

% the aibrake_sim.slx should be run once beforehand in order for this
% script to work without issues

% bounds of the decision variables can be selected by modifiying VarMin and
% VarMax, while iteration number can be modified with MaxIt

% the variables output should be saved to the pid_var.mat file manually
% using the save command commented at the end of this script only once the
% tuning process is complete

% further information on how to tune the PID parameters can be found in the
% README.md located inside the Airbrakes folder

%% Problem Definiton

%%% PLEASE CHOOSE APPROPRIATE BOUNDS%%%

nVar = 4;        % Number of Unknown (Decision) Variables
VarSize = [1 nVar];         % Matrix Size of Decision Variables
VarMin = [18 0.2 9.5 50];	% Lower Bound of Decision Variables
VarMax = [19 0.35 10.5 54];    % Upper Bound of Decision Variables
%% Parameters of PSO
MaxIt = 10;     % Maximum Number of Iterations
nPop = 4;      % Population Size (Swarm Size)
w = 1;          % Intertia Coefficient
wdamp = 0.99;   % Damping Ratio of Inertia Coefficient
c1 = 2;         % Personal Acceleration Coefficient
c2 = 2;         % Social Acceleration Coefficient
% The Flag for Showing Iteration Information
MaxVelocity = 0.2*(VarMax-VarMin);
MinVelocity = -MaxVelocity;
%% Initialization
% The Particle Template
empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];
% Create Population Array
particle = repmat(empty_particle, nPop, 1);
% Initialize Global Best
GlobalBest.Cost = inf;
% Initialize Population Members
for i=1:nPop
    % Generate Random Solution
    particle(i).Position = unifrnd(VarMin, VarMax, VarSize);
    P = particle(i).Position(1);
    I = particle(i).Position(2);
    D = particle(i).Position(3);
    N = particle(i).Position(4);
    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    % Evaluation
    out = sim('Airbrake_sim');
    particle(i).Cost = out.Fobj(end,1);
    % Update the Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    % Update Global Best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
end

% Array to Hold Best Cost Value on Each Iteration
BestCosts = zeros(MaxIt, 1);
%% Main Loop of PSO
for it=1:MaxIt
    for i=1:nPop
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position) ...
            + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity, MinVelocity);
        particle(i).Velocity = min(particle(i).Velocity, MaxVelocity);
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        % Apply Lower and Upper Bound Limits
        particle(i).Position = max(particle(i).Position, VarMin);
        particle(i).Position = min(particle(i).Position, VarMax);
        P = particle(i).Position(1);
        I = particle(i).Position(2);
        D = particle(i).Position(3);
        N = particle(i).Position(4);
        % Evaluation
        out = sim('Airbrake_sim');
        particle(i).Cost = out.Fobj(end,1);
        % Update Personal Best
        if particle(i).Cost < particle(i).Best.Cost
            
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
        end
    end
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCosts(it))]);
    % Damping Inertia Coefficient
    w = w * wdamp;
end
P = GlobalBest.Position(1)
I = GlobalBest.Position(2)
D = GlobalBest.Position(3)
N = GlobalBest.Position(4)


open('Airbrake_sim');
sim('Airbrake_sim');

% save pid_var.mat P I D N