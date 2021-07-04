function [classAdj_means, classAdj_stdevs, totalTrees_mean, totalTrees_stdev] = OS_MultivariateRuns_varianceByPermuting(MST, classOccupancy, iter)

% 1. set up some variables
N = sum(classOccupancy);
k = length(classOccupancy);
classAdjacencies = zeros(k, k, iter);
totalTrees = zeros(iter, 1);
classLabelsInit = []; %                make class labels vector
for u = 1:k
    classLabelsInit = [classLabelsInit; u*ones(classOccupancy(u), 1)]; 
end

% 2. main loop: permute labels, get MST class adjacencies with new labels, repeat
parfor q = 1:iter
    % 2.1 random permutation of the labels
    classLabels = classLabelsInit(randperm(N)); 
    
    % 2.2 loop through MST.Edges.EndNodes, adding 1 to each corresponding edge type in classAdj
    classAdj = zeros(k, k);
    edgeClasses = [classLabels(MST.Edges.EndNodes(:, 1)), classLabels(MST.Edges.EndNodes(:, 2))];
    for v = 1:N-1
        classAdj(edgeClasses(v, 1), edgeClasses(v, 2)) = classAdj(edgeClasses(v, 1), edgeClasses(v, 2)) + 1;
    end
    
    classAdjacencies(:, :, q) = classAdj + classAdj' - diag(diag(classAdj)); %     symmetrize and add to the stack
    totalTrees(q) = N - trace(classAdj); %                     get total trees produced, by subtracting trace from N
    
end % 
fprintf('\n')

% 3. get class adjacency mean and stdev
classAdj_means = mean(classAdjacencies, 3);
classAdj_stdevs = std(classAdjacencies, 0, 3);
totalTrees_mean = mean(totalTrees);
totalTrees_stdev = std(totalTrees);

end