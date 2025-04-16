function kFunction = computeRipleysK(points, rRange)
    numPoints = size(points, 1);
    kFunction = zeros(size(rRange));
    
    for i = 1:length(rRange)
        r = rRange(i);
        
        % Count the number of points within distance r for each point
        numNeighbors = zeros(numPoints, 1);
        for j = 1:numPoints
            distances = sqrt(sum((points - points(j, :)).^2, 2));
            numNeighbors(j) = sum(distances <= r);
        end
        
        % Compute the average number of neighbors
        avgNumNeighbors = mean(numNeighbors);
        
        % Compute the K function value
        kFunction(i) = avgNumNeighbors / numPoints;
    end
end