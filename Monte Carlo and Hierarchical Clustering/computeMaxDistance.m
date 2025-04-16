function maxDistance = computeMaxDistance(points)
    numPoints = size(points, 1);
    maxDistance = 0;
    
    for i = 1:numPoints
        for j = i+1:numPoints
            % Compute the Euclidean distance between points i and j
            distance = norm(points(i, :) - points(j, :));
            
            % Update maxDistance if the current distance is greater
            if distance > maxDistance
                maxDistance = distance;
            end
        end
    end
end