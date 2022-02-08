function infectedRMSE = RMSEFinder(alpha0)
global A seedNodes 
% Import the network and get a bit of data about it
% A = importdata("highSchool.mat");
% seedNode = 25;
numNodes = size(A,1);
[edgeRows, edgeCols] = find(A);
edgeArray = [edgeRows edgeCols];

% I have to manually input these things. If I pass them through the
% function above then fminsearch will try to optimise the lambdas. I'm sure
% there's a cleverer way of doing this but I'm too lazy to figure it out.
lambdai = 1;
lambdar = 0.1;
params = [lambdai,lambdar,alpha0];

T_max = 200;
timeStep = 0.05*lambdai;

tspan = [0:timeStep:T_max];

% I have always initially infected node 1, so initial conditions hard coded
% once again
initConds = [ones(numNodes,1),zeros(numNodes,1)];
initConds(seedNodes,1) = 0;
initConds(seedNodes,2) = 1;
y0 = [initConds(:,1); initConds(:,2)];

%Just a check in case something screwy happens
if alpha0 >1
    disp("Error, alpha > 1")
end
if alpha0 < 0
    disp("Error, alpha < 0")
end
%[t,pHybrid] = ode23tb(@(t,pAlpha) hybrid(t,pAlpha,numNodes,A,alpha),tspan,y0);

% Call hybrid method
[t,yHybrid] = ode45(@(tForFn,y) ...
    sirSpecDE(tForFn,y,params,edgeArray(:,1),edgeArray(:,2)), ...
    tspan, y0);

%% Get averages from Hybrid
sSol = yHybrid(:,1:numNodes);
iSol = yHybrid(:,numNodes+1:2*numNodes);
rSol = 1-iSol-sSol;

% averagesSol = sum(sSOl')/numNodes;
averageISol = sum(iSol')/numNodes;
averageRSol = sum(rSol')/numNodes;

%% Import averages from Numerical simulations
avgGilI = importdata("numResults.mat");

%I print this just to keep track
alpha0 

% Calculating the error:
infectedRMSE= sqrt(sum((avgGilI-averageISol).^2));

