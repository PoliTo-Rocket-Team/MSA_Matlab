clc;
% clear;
close all;

% WARNING
% This script requires atleast a full run of the airbrake_sim.slx simulator
% to work accordingly

% INFO
% this PSO script is used for generating the P, I, D and N values used in
% the PID which controls the air brakes control system

% the aibrake_sim.slx needs to be run in order for this script to work

% bounds of the decision variables can be selected by modifiying VarMin and
% VarMax

% the variables output must be transcribed manually into the PID block
% inside airbrake_sim.slx (but they may be temporarely changed to the
% output variables themselves to speed up the optimization process)

%% Problem Definiton

%%% PLEASE CHOOSE APPROPRIATE BOUNDS%%%

nVar = 4;        % Number of Unknown (Decision) Variables
VarSize = [1 nVar];         % Matrix Size of Decision Variables
VarMin = [1 1 1 1];	% Lower Bound of Decision Variables
VarMax = [100 10 10 100];    % Upper Bound of Decision Variables
%% Parameters of PSO
MaxIt = 3;     % Maximum Number of Iterations
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
    sim('airbrake_sim');
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
        sim('airbrake_sim');
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
open('airbrake_sim');
sim('airbrake_sim');