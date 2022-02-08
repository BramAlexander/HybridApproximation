%% sirHybridModel.m
% 2020-04-08

function [...
    sSol, ...               Solution for susceptibles (numNodes by numTimes)
    iSol, ...               Solution for infecteds (numNodes by numTimes)
    rSol ...                Solution for recovereds (numNodes by numTimes)
    ] = sirHybridModel(...
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    params, ...             Model parameters [lambda, gamma, alpha]
    initConds, ...          Initial conditions for all nodes (numNodes by 2)
    t ...                   Vector of times (1 by numTimes)
    )
    
% Notes on input arguments

% edgeArray
% Each row of edgeArray corresponds to an edge of the graph where the node
% indicated by the first column can spread disease to the node indicated by
% the second column. For an undirected graph, there are two rows of
% edgeArray for every edge. edgeArray can be generated from the adjacency
% matrix using find.

% initConds
% Each row of initConds gives the probability (between 0 and 1) of a node
% being susceptible (first column) or infected (second column) at time 0.


%% Extracting additional parameters

% Number of nodes in system
numNodes = size(initConds,1);

% Initial conditions
y0 = [initConds(:,1); initConds(:,2)];

% Maximum time
maxTime = t(end);

%% Solve differential equations

% Solve sirSpecDE using ODE45
sol = ode45(@(tForFn,y) ...
    sirSpecDE(tForFn,y,params,edgeArray(:,1),edgeArray(:,2)), ...
    [0 maxTime], y0);

% Evaluate solution at specified times
ySol = deval(sol,t);

% Convert solution to solutions for susceptibles, infecteds, recovereds
sSol = ySol(1:numNodes,:);
iSol = ySol((numNodes+1):2*numNodes,:);
rSol = 1 - sSol - iSol;

end

%% Nested function to specify differential equation.
function dydt = sirSpecDE(~,y,params,mainNode,neighbourNode)

% Note:
% (mainNode,neighbourNode) corresponds to the (i,j) coefficient of one of
% the 1s in the adjacency matrix. They are used in treeModelRates and
% classicModelRates in order to avoid doing too many unnecessary
% calculations of things that will be multiplied by zero.

%% Preliminaries

% Parameter extraction
numNodes = numel(y)/2;
lambda = params(1);
gamma = params(2);
alpha = params(3);

% Initialisation
dydt = zeros(size(y));
treeModelRates = zeros(numNodes);
classicModelRates = zeros(numNodes);
s = y(1:numNodes);
i = y((numNodes+1):2*numNodes);

%% Calculate rates

% In order to calculate the rate of infection at each node, we construct
% matrices (initialised above) where the (i,j) element corresponds to the
% rate of infection at i due to a neighbour at j. These matrices are only
% nonzero where there is an edge between i and j, and so the calculations
% are only done where [i,j] = [mainNode,neighbourNode], based on the
% nonzero elements of the adjacency matrix.

% adjIndex is the index (in single number form) of the nonzero elements of
% the adjancency matrix.
adjIndex = (neighbourNode-1)*numNodes+mainNode;

% Infection rate according to the model that is exact for trees.
treeModelRates(adjIndex) = ...
    max((lambda+gamma)*s(mainNode) - lambda*s(neighbourNode) - gamma,0);

% Infection rate based on treating nodes as independent.
classicModelRates(adjIndex) = ...
    lambda*s(mainNode).*i(neighbourNode);

% Infection rate from a weighted average of the two models.
totalInfRate = (1-alpha)*sum(treeModelRates,2) ...
    + alpha*sum(classicModelRates,2);

% Rate of change of susceptible probability is -totalInfRate.
% Rate of change of infected probability is +totalInfRate - recoveryRate.
dydt(1:numNodes) = -totalInfRate;
dydt((numNodes+1):2*numNodes) = totalInfRate - gamma*i;

end