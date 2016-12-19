function [finalPopulation] = GeneticSimulation
% ==============================================================
% PROTOTYPE SIMULATION FOR GAME THEORY PROJECT
% in the course
% Game Theory and Rationality ENM140, Chalmers
% =============================================================
%
% The simulation features an evoulutionary all-vs-all scenario 
% with a population of players playing 'The Traffic Game' for a
% set number of paths and cost parameter.
%
% The 'learning algorithm' chosen is similar to the type of
% evolution found in a standard Genetic Algorithm for
% optimization.
%
% By Simon Nilsson (simnilss)
% Last updated 2016-12-08


% =========== GAME PARAMETERS =============================
populationSize = 50;    % N in mathematical description
nActions = 20;          % m in mathematical description
costParameter = 1;      % c in mathematical description
nRoundsPerGeneration = 100;

% ======== SIMULATION PARAMETERS ==========================
nGenerations = 1000;
crossoverProb = 0.1;
mutationProb = 1/(populationSize*nActions);
criticalMutationProb = 0.1;
mutationRate = 0.05;
% ---------------------------------------------------------

% Set seed based on current time for RNG
rng('shuffle');

% ======= VARIABLES AND 'ALLOCATION' ======================
% INITIAL STRATEGIES
% Generate random initial population strategies
population = NormalizeProbabilities(rand(populationSize, nActions));

% Figures and plot handles
figure(1); clf;
stratPlot = [];
for i=1:min(4, populationSize)
    subplot(2,2,i);
    stratPlot(i) = bar(population(i,:));
    ylim([0 1])
    title(['Strategy ' num2str(i)])
    xlabel('Path')
    ylabel('Probability')
end

figure(2); clf;
meanScore = zeros(nGenerations,1);
scorePlot = plot(meanScore);
title('Mean score')
xlabel('Generation')
ylabel('Score')

% ============= EVOLUTIONARY SIMULATION ==============
% ----------------------------------------------------

bestIndividual = zeros(1, nActions);
for iGeneration = 1:nGenerations
    % ============ PLAY THE GAMES ===============
    score = zeros(1, populationSize);
    for iRound = 1:nRoundsPerGeneration
        % Generate random numbers for path choosing
        r = rand(populationSize, 1);
        % Compute random number thresholds
        thresholds = cumsum(population, 2);


        % Find what paths are chosen
        % -----------------------------
        % This sorts the logical array in descending order and
        % storing the index order for each row, i.e. the index of
        % the first occurrence of a 1 will be in the first column
        % of paths
        % Player i's path choice will be at index i of 'paths'
        [~, pathChoice] = sort(thresholds > repmat(r, [1,size(thresholds,2)]), 2, 'descend');
        pathChoice = pathChoice(:,1);


        % Count occurrences of path (path(i) in [1,nActions]
        playersOnPath = zeros(nActions, 1);
        for i=1:nActions
            playersOnPath(i) = sum(pathChoice == i);
        end
    
        % ---- Evaluate scores -----
        for i = 1:populationSize
           score(i) = score(i) + nActions+1-pathChoice(i) ...
               - costParameter*(playersOnPath(pathChoice(i)) - 1);
        end
    end
    % Compute mean score for this generation
    meanScore(iGeneration) = mean(score) / nRoundsPerGeneration;
    
    % Compute the best result and which player got it
    [~, bestResultIndex] = max(score);
    bestIndividual = population(bestResultIndex,:);
    
    % Temporarily assign next generation as the current one
    nextGeneration = population;
    
    % ========== SELECTION AND CROSSOVER ===========
    for i=1:2:populationSize
        
        % Pick two individuals
        individual_1 = Selection(population, score);
        individual_2 = Selection(population, score);
        
        % Check if crossover occurs
        if (rand < crossoverProb)
            [individual_1, individual_2] = Crossover(individual_1, individual_2);
        end
       
        % Add picked individuals to next generation
        nextGeneration(i,:) = individual_1;
        
        % Ugly, but without this the population increases by one for odd
        % population sizes
        if (i < populationSize)
            nextGeneration(i+1,:) = individual_2;
        end
    end
    
    % ========= MUTATION =========
    nextGeneration = Mutate(nextGeneration, mutationProb, ...
        mutationRate, criticalMutationProb);
    
    % Add best individual to population so that it doesn't disappear
    nextGeneration(1,:) = bestIndividual;
    % Switch to the new population
    population = nextGeneration;
    
    
    % ========== UPDATE PLOTS ===========
    % Update strategy plot
    for i = 1:min(4, populationSize)
        set(stratPlot(i), 'YData', population(i,:));
    end
    % Update mean score
    set(scorePlot, 'YData', meanScore);
    
    % Arbitrary pause for plots to update
    pause(0.05)
    
end

% Return final population
finalPopulation = population;

end

% ============ HELPER FUNCITON DEFINITIONS ===============
% --------------------------------------------------------

function probs = NormalizeProbabilities(P)
    
    s = sum(P, 2);
    probs = P ./ repmat(s, [1, size(P,2)]);
end

function [outGene1, outGene2] = Crossover(gene1, gene2)

    minLength = min(length(gene1), length(gene2)) - 1;
    splitPoint = randi(minLength);
    outGene1 = [gene1(1:splitPoint) gene2(splitPoint+1:end)];
    outGene2 = [gene2(1:splitPoint) gene1(splitPoint+1:end)];
end

function outPopulation = Mutate(population, mutProb, rate, criticalMutProb)
    
    % Generate mutations in the genome
    mutations = rand(size(population)) < mutProb;
    % Evaluate how many of them are 'critical mutations'
    criticalMutations = logical((rand(size(mutations)) < criticalMutProb) .* mutations);
    
    % Carry out mutations and critical mutations
    outPopulation = population;
    outPopulation(mutations) = outPopulation(mutations) ...
        + rate * (2*rand(sum(sum(mutations)), 1) - 1);
    outPopulation(criticalMutations) = rand(sum(sum(criticalMutations)), 1);
    
    % Make sure there are no negative probabilities
    negative = outPopulation < 0;
    outPopulation(negative) = 0;
    
    % Renormalize the probabilities
    outPopulation = NormalizeProbabilities(outPopulation);
    
end

function individual = Selection(population, score)
    
    % Compute score thresholds for roulette wheel selection
    scoreThreshold = cumsum(score) / sum(score);
    
    % Use Roulette-wheel selection
    r = rand;    
    temp = repmat(r, [1, size(scoreThreshold, 2)]) ...
        < scoreThreshold;
    choice = find(temp, 1);
    
    individual = population(choice,:);
    
end
