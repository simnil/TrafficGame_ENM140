function [] = SpatialSimulation
% ===================================================================
% SPATIAL PROTOTYPE SIMULATION FOR GAME THEORY PROJECT
% Traffic Coordination Game Project, in the course
% Game Theory and Rationality ENM140, Chalmers
% ===================================================================
%
% This simulation add a spatial aspect to the Traffic Coordination 
% Game studied. The agents on the grid play only against their
% neighbours, i.e. participate in 9 games (including their own) when
% using the von-Neumann neighbourhood. Optionally the Moore
% neighbourhood can be used.
%
% Only pure strategies are used and the update model is that agents
% choose the strategy (path) of the agent in the neighbourhood with
% the highest score. They will prioritize their own strategy in this.
%
% By Simon Nilsson (simnilss)
% Last update 2016-12-12



% ===== SIMULATION PARAMETERS ===========================================
gridSize = 200; % The number of agents will be equal to gridSize*gridSize
nPaths = 3;     % The number of paths in each subgame (1..neighbours-1)
cost = 1;       % The c parameter of the mathematical description
nTimesteps = 2e2;
% The probability that the agent switches path regardless of score
% (to introduce perturbations)
switchProb = 10/(gridSize*gridSize);

% ===== INITIALIZE VARIABLES ============================================
% Initialize a population of gridSize*gridSize agents on a grid with
% random pure strategies represented by integers in the range [1, nPaths]
%population = randi(nPaths, [gridSize*gridSize, 1]);

% OPTION: initialize the entire population with the same strategy
% and let either a single player start with a different one or let
% noise perturbations initiate evolution.
population = ones(gridSize*gridSize, 1);
%population(floor(3*gridSize*gridSize/4)) = 2;

meanScore = zeros(nTimesteps, 1);
customColors =  [0 0 0.5;
                0 0 1;
                1 1 0.1;
                0 0.75 0.125;
                0.6 0 0.6;
                0.8 0.1 0.5;
                1 0 0];

% -- Create VideoWriter struct and edit settings ---
video = VideoWriter('spatial result3.avi');
video.Quality = 100;
open(video);


% ---- Figures ----------------
figure(1); clf;
plotHandle = image( reshape(population,[gridSize,gridSize]) );
% Custom colormap
colormap(customColors(1:nPaths,:));
% Black and white colormap
%colormap('gray')
colorbar

figure(2); clf;
scorePlot = plot(1:nTimesteps, meanScore);
xlabel('Timestep')
ylabel('Mean score')

% ============= SIMULATION LOOP =============
% -------------------------------------------
for t = 1:nTimesteps
    
    % Let everyone participate in their neighbourhood games
    score = zeros(size(population));
    for i = 1:length(population)
       score = score + SubGame(population, i, gridSize, nPaths, cost);
    end
    meanScore(t) = mean(score);
    
    % Update everyone's strategy
    newStrategies = zeros(size(population));
    for i = 1:length(population)
        newStrategies(i) = UpdateStrategy(population, score, i, gridSize);
    end
    switches = rand(size(newStrategies)) < switchProb;
    newStrategies(switches) = randi(nPaths, size(newStrategies(switches)));
    population = newStrategies;
    
    
    % Update plots
    % ------------
    % NOTE!: But update only every other timestep, for a more 'nice-looking'
    % time evolution. Rapid fluctuations cause visual 'stutter' which is
    % unpleasant to watch.
    if (mod(t, 2) == 0)
        img = reshape(population, [gridSize, gridSize]);
        set(plotHandle, 'CData', img);
        % arbitrary pause for graphics to update
        pause(0.05);
        frame = im2frame(img, customColors);
        writeVideo(video, frame);
        set(scorePlot, 'YData', meanScore);
        % arbitrary pause for graphics to update
        pause(0.05);
    end
end

close(video);

end

function score = SubGame(population, coord, gridSize, nPaths, cost)
    
    % Find neighbours, with wrap-around
    neighbours = Neighbourhood(coord, gridSize);
    
    % Extract player strategies and compute path distribution
    strategies = population(neighbours);
    playersOnPath = histcounts(strategies, 0.5:nPaths+0.5)';
    
    score = zeros(size(population));
    score(neighbours) = nPaths+1 - strategies - cost*(playersOnPath(strategies)-1);
    
end


% Returns a vector containing the indices corresponding to the neighbouring
% tiles in the grid (with wrap-around)
function neighbours = Neighbourhood(coord, gridSize)

    left        = mod( coord-2, gridSize*gridSize ) + 1;
    right       = mod( coord, gridSize*gridSize ) + 1;
    up          = mod( coord-1 - gridSize, gridSize*gridSize ) + 1;
    down        = mod( coord-1 + gridSize, gridSize*gridSize ) + 1;
    upleft      = mod( coord-2 - gridSize, gridSize*gridSize ) + 1;
    upright     = mod( coord - gridSize, gridSize*gridSize ) + 1;
    downleft    = mod( coord-2 + gridSize, gridSize*gridSize ) + 1;
    downright   = mod( coord + gridSize, gridSize*gridSize ) + 1;
    
    % von-Neumann neighbourhood
    neighbours = [coord left upleft up upright right downright down downleft]';
    % OPTION:
    % Moore neighbourhood
    %neighbours = [coord left up right down]';
end

function strategy = UpdateStrategy(population, score, coord, gridSize)

    neighbours = Neighbourhood(coord, gridSize);
    
    % Find who got the best scores
    [scores, index] = sort(score(neighbours), 1, 'descend');
    
    % If player was one of them, keep the current strategy
    if (score(coord) == scores(1))
        strategy = population(coord);
        return;
    
    % If someone else was better, choose one of those strategies
    else
        % How many got the best score
        nBest = sum(scores == scores(1));
        % Choose one of them at random
        choice = neighbours(index(randi(nBest)));
        strategy = population( choice );
        return;
    end

end