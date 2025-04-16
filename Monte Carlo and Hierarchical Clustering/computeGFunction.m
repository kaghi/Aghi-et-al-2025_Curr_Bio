function gFunction = computeGFunction(points, rRange)
    numPoints = size(points, 1);
    gFunction = zeros(size(rRange));
    
    for i = 1:length(rRange)
        r = rRange(i);
        
        % Compute the number of neighbors within distance r for each point
        numNeighbors = zeros(numPoints, 1);
        for j = 1:numPoints
            distances = sqrt(sum((points - points(j, :)).^2, 2));
            numNeighbors(j) = sum(distances > 0 & distances <= r);
        end
        
        % Compute the average number of neighbors
        avgNumNeighbors = mean(numNeighbors);
        
        % Compute the G function value
        gFunction(i) = avgNumNeighbors / numPoints;
    end
end
