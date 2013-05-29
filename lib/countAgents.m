function grid = countAgents (agents,nCells) 
% COUNTAGENTS count the number the of agents per cell of the grid

grid = zeros(nCells);

classes =linspace(0,1,nCells);


for i=1:size(agents,2)
    % Columns
    idxC = find((classes <= agents(1,i)), 1, 'last');
    % Rows 
    idxR = nCells+1 - find((classes <= agents(2,i)), 1, 'last');
    grid(idxR,idxC) = grid(idxR,idxC) + 1; 
end

