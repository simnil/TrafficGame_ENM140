function [] = Schelling
% ==============================================================
% SCHELLING MODEL SIMULATION
% in the course
% Game Theory and Rationality ENM140, Chalmers
% ==============================================================
% By Simon Nilsson (simnilss)
% Last updated 2017-01-01

% =================== SIMULATION PARAMETERS ====================
gridSize = 128;
satisfyThreshold = 0.7;
emptyFrac = 0.1;
nTimesteps = 5e2;

customColors = [1 1 1;
                1 0 0;
                0 0 1];


% --- Initialize agents ---
population = ones(gridSize*gridSize, 1);
population( rand(size(population)) < 0.5 ) = 2;
% Make some tiles empty, so movement is possible
population( rand(size(population)) < emptyFrac) = 0;
emptyIdx = find(population == 0);

% Pre-allocate satisfaction vector
satisfaction = zeros(size(population));


% Set up plot
figure(1); clf;
imgHandle = imagesc( reshape(population, [gridSize,gridSize]) );
colormap(customColors);
colorbar


% ============= SIMULATION LOOP =============
% -------------------------------------------
for iTimestep=1:nTimesteps
    
    % -- Update the satisfaction of all agents --
    satisfaction(population == 0) = 1;
    for iAgent=1:gridSize*gridSize
        if ( population(iAgent) )
            satisfaction(iAgent) = CheckSatisfaction(iAgent, satisfyThreshold, ...
                population, gridSize);
        end
    end
    
    % -- Move unsatisfied agents --
    [population, emptyIdx] = MoveUnsatisfied(population, satisfaction, emptyIdx);
    
    % -- Update graphics in plot --
    img = reshape(population, [gridSize,gridSize]);
    set(imgHandle, 'CData', img);
    % Arbitrary pause for graphics to update
    pause(0.01)
    
    % Break if no-one is unsatisfied, as nothing will happen
    if (sum(satisfaction == 0) == 0)
        break
    end
end

end

% ================ HELPER FUNCTIONS DEFINITIONS ================
% --------------------------------------------------------------

% Returns a vector containing the indices corresponding to the neighbouring
% tiles in the grid (with wrap-around)
function neighbours = Neighbourhood(idx, gridSize)

    left        = mod( idx-2, gridSize*gridSize ) + 1;
    right       = mod( idx, gridSize*gridSize ) + 1;
    up          = mod( idx-1 - gridSize, gridSize*gridSize ) + 1;
    down        = mod( idx-1 + gridSize, gridSize*gridSize ) + 1;
    upleft      = mod( idx-2 - gridSize, gridSize*gridSize ) + 1;
    upright     = mod( idx - gridSize, gridSize*gridSize ) + 1;
    downleft    = mod( idx-2 + gridSize, gridSize*gridSize ) + 1;
    downright   = mod( idx + gridSize, gridSize*gridSize ) + 1;
    
    % von-Neumann neighbourhood
    neighbours = [left upleft up upright right downright down downleft]';
    % OPTION:
    % Moore neighbourhood
    %neighbours = [left up right down]';
end

function satisfaction = CheckSatisfaction(idx, threshold, population, gridSize)
    
    agent = population(idx);
    neighbours = population( Neighbourhood(idx, gridSize) );
    % Check how many have the same type
    nSame = sum( neighbours == agent );
    % Obtain fraction of neighbours with same type (ignoring empty tiles)
    fracSame = nSame / sum( neighbours ~= 0 );
    
    satisfaction = fracSame >= threshold;
end

function [newPopulation, newEmptyIdx] = MoveUnsatisfied(population, satisfaction, emptyIdx)

    newPopulation = population;
    
    % Find unsatisfied agents
    unsatisfied = find(satisfaction == 0);
    nUnsatisfied = length(unsatisfied);
    % Move them one by one
    for i = 1:nUnsatisfied
        
        % Agent's current location
        agent = unsatisfied(i);
        
        type = newPopulation(agent);
        newPopulation(agent) = 0;
        % Find new location
        temp = randi(length(emptyIdx));
        newLoc = emptyIdx(temp);
        % Remove index from list of empty ones
        emptyIdx(temp) = [];
        
        % Move agent to the new location
        newPopulation(newLoc) = type;
        % Add old location to list of empty ones
        emptyIdx = [emptyIdx; agent];
        
    end
    
    newEmptyIdx = emptyIdx;
    
end