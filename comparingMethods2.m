clear all;
close all;
%%
global A seedNodes
seedNodes = 25;
plt = [1,7,8,27,34];
plt = plt(plt~=seedNodes);
%%
A = importdata("karate.mat");
[edgeRows, edgeCols] = find(A);
edgeArray = [edgeRows edgeCols];
%Initial Conditions
numNodes = size(A,2);

%seedNodes = randi(numNodes,1,floor(0.05*numNodes));
% y0 = ones(2,numNodes);
% y0(2,:) = 0;
% y0(2,seedNodes) = 1;

initConds = [ones(numNodes,1),zeros(numNodes,1)];
initConds(seedNodes',1) = 0;
initConds(seedNodes',2) = 1;

%Setting up time

T_max = 200;
timeStep = 0.05;

%Parameters

lambdaI = 1;
lambdaR = .1;
alpha = 0.5;

params = [lambdaI,lambdaR,alpha];
paramsMFA = [lambdaI,lambdaR,1];
paramsMin = [lambdaI,lambdaR,0];

tspan = [0:timeStep:T_max];

k = sum(A);
k = k/(max(k));
% Calling methods
%% Call min(AB)
% pMFA = MFA_SI(A,numNodes,T_max,y0,lambdai);
numRuns = 10000;
tic
[probS,probI,probR] = sirGillespie(A,params(1:2),initConds,tspan,numRuns);
toc

tic
[Time, pairResults] = pair_based_model(A,lambdaR,seedNodes,tspan);
%pairResults = pairApprox_SIR(A,T_max,lambdaI,lambdaR,timeStep);
toc

tic
[sSol,iSol,rSol] = ...
    sirHybridModel(...
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    params, ...             Model parameters [lambda, gamma, alpha]
    initConds, ...          Initial conditions for all nodes (numNodes by 2)
    tspan ...                   Vector of times (1 by numTimes)
    );
toc
% [sSolK,iSolK,rSolK] = ...
%     sirHybridModelK(...
%     edgeArray, ...          List of edges (numDirectedEdges by 2)
%     params, ...             Model parameters [lambda, gamma, alpha]
%     initConds, ...          Initial conditions for all nodes (numNodes by 2)
%     tspan, ...                   Vector of times (1 by numTimes)
%     k ...
%     );
tic
[sMFA,iMFA,rMFA] = ...
    sirHybridModel(...
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    paramsMFA, ...             Model parameters [lambda, gamma, alpha]
    initConds, ...          Initial conditions for all nodes (numNodes by 2)
    tspan ...                   Vector of times (1 by numTimes)
    );
toc

tic
[sMin,iMin,rMin] = ...
    sirHybridModel(...
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    paramsMin, ...             Model parameters [lambda, gamma, alpha]
    initConds, ...          Initial conditions for all nodes (numNodes by 2)
    tspan ...                   Vector of times (1 by numTimes)
    );
toc
avgGilS = sum(probS)/numNodes;
avgGilI = sum(probI)/numNodes;
avgGilR = sum(probR)/numNodes;

save('numResults.mat','avgGilI');

avgPairS = sum(pairResults(:,1:numNodes)')/numNodes;
avgPairI = sum(pairResults(:,numNodes+1:2*numNodes)')/numNodes;
avgPairR = 1 - avgPairS - avgPairI;
%avgPairR = sum(pairResults(:,201:300)')/numNodes;
 
avgHybridS = sum(sSol)/numNodes;
avgHybridI = sum(iSol)/numNodes;
avgHybridR = sum(rSol)/numNodes;

% avgHybridSK = sum(sSolK)/numNodes;
% avgHybridIK = sum(iSolK)/numNodes;
% avgHybridRK = sum(rSolK)/numNodes;

avgMFAS = sum(sMFA)/numNodes;
avgMFAI = sum(iMFA)/numNodes;
avgMFAR = sum(rMFA)/numNodes;

avgMinS = sum(sMin)/numNodes;
avgMinI = sum(iMin)/numNodes;
avgMinR = sum(rMin)/numNodes;

optimAlpha = fminsearch(@RMSEFinder,0.5);
optimParams = [lambdaI,lambdaR,optimAlpha];

tic
[sSolOptim,iSolOptim,rSolOptim] = ...
    sirHybridModel(...
    edgeArray, ...          List of edges (numDirectedEdges by 2)
    optimParams, ...             Model parameters [lambda, gamma, alpha]
    initConds, ...          Initial conditions for all nodes (numNodes by 2)
    tspan ...                   Vector of times (1 by numTimes)
    );
toc

avgOptimHybridS = sum(sSolOptim)/numNodes;
avgOptimHybridI = sum(iSolOptim)/numNodes;
avgOptimHybridR = sum(rSolOptim)/numNodes;

figure
plot(tspan,avgGilS,'o')
hold on
plot(tspan,avgHybridS)
plot(tspan,avgOptimHybridS)
plot(tspan,avgPairS,'--')
plot(tspan,avgMFAS,'-x')
plot(tspan,avgMinS,'-.')
hold off
legend("Numerical","Hybrid","OptimaAlpha","Pair","MFA","Tree")
title("Susceptible")

figure
plot(tspan,avgGilI,'o')
hold on
plot(tspan,avgHybridI,'+')
plot(tspan,avgOptimHybridI,'k')
plot(tspan,avgPairI,'--')
plot(tspan,avgMFAI,'-x')
plot(tspan,avgMinI,'-.')
hold off
legend("Numerical","Hybrid","OptimaAlpha","Pair","MFA","Tree")
ylabel("Proportion Infected")
xlabel("Time")
figure
plot(tspan,avgGilR,'o')
hold on
plot(tspan,avgHybridR)
plot(tspan,avgOptimHybridR)
plot(tspan,avgPairR,'--')
plot(tspan,avgMFAR,'-x')
plot(tspan,avgMinR,'-.')
legend("Numerical","Hybrid","OptimaAlpha","Pair","MFA","Tree")
title("Recovered")

% figure
% plot(tspan,1-avgGilS,'o')
% hold on
% plot(tspan,1- avgPairS,'--')
% plot(tspan,1 - avgHybridS)
% plot(tspan,1-avgMFAS,'-x')
% plot(tspan,1-avgMinS,'-.')
% hold off
% legend("Numerical","Pair","Hybrid","MFA","Tree")
% title("1 - Susceptible")

infectedMFAError = sqrt(sum((avgGilI-avgMFAI).^2))
recoveredMFADifference = avgGilR(end)-avgMFAR(end);
notSusceptibleMFAError = sqrt(sum((1-avgGilS)-(1-avgMFAS)).^2);

infectedMinError = sqrt(sum((avgGilI-avgMinI).^2));
recoveredMinDifference = avgGilR(end)-avgMinR(end);
notSusceptibleMinError = sqrt(sum((1-avgGilS)-(1-avgMinS)).^2);

alpha
infectedHybridError = sqrt(sum((avgGilI-avgHybridI).^2))
recoveredHybridDifference = avgGilR(end)-avgHybridR(end);
notSusceptibleHybridError = sqrt(sum((1-avgGilS)-(1-avgHybridS)).^2);

optimAlpha
infectedHybridOptimError = sqrt(sum((avgGilI-avgOptimHybridI).^2))
recoveredHybridOptimDifference = avgGilR(end)-avgOptimHybridR(end);
notSusceptibleHybridOptimError = sqrt(sum((1-avgGilS)-(1-avgOptimHybridS)).^2);

infectedPairError = sqrt(sum((avgGilI-avgPairI).^2))
recoveredPairDifference = avgGilR(end)-avgPairR(end);
notSusceptiblePairError = sqrt(sum((1-avgGilS)-(1-avgPairS)).^2);

[MaxMFAI,IMFA] = max(avgMFAI);
[MaxMinI,IMin] = max(avgMinI);
[MaxPairI,IPair] = max(avgPairI);
[MaxHybridI,IHybrid] = max(avgHybridI);
[MaxGilI,IGil] = max(avgGilI);

MFAPeakError = MaxMFAI - MaxGilI;
MFATImeError = IMFA - IGil;

MinPeakError = MaxMinI - MaxGilI;
MinTImeError = IMin - IGil;

PairPeakError = MaxPairI - MaxGilI;
PairTImeError = IPair - IGil;

HybridPeakError = MaxHybridI - MaxGilI;
HybridTImeError = IHybrid - IGil;


figure
plot(probI(plt,:)','o')
hold on
plot(iSol(plt,:)')
legend([string(plt),string(plt)])

figure
plot(probS(plt,:)','o')
hold on
plot(sSol(plt,:)')
legend([string(plt),string(plt)])
% RMSEPair = sqrt(sum(sum(((gillespieA(:,:,1)-pairResults(:,1:numNodes).^2)*lambdaI))));
%  RMSEHybrid = sqrt(sum(sum(((gillespieA(:,:,1)-sSol').^2)*lambdaI)))
% % 
