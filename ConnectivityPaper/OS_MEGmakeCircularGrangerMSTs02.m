
%% 1. Set some variables
tic
clc
clear
band = 'delta';
dt = 1; %                                       Time offset (integer, 0-2; # of 1.5-second trials)
symmetrizeMethod = 'max'; %         Distance matrix symmetrization: 'max' takes the larger abs(correlation) at positions ij and ji; 'avg' averages the two 
subtractParcAdj = true; %                For dt = 0, if 'true', we subtract out the parcel adjacency effect

indicesNMP = {[1:35], [36:72], [73:89]};
labelvecNMP = [ones(1, length(indicesNMP{1})), 2*ones(1, length(indicesNMP{2})), 3*ones(1, length(indicesNMP{3}))];

classes = {'Normal', 'MCI', 'PDD'};

%% 2. Get coordinates for 116 points in a circle, and display with parcel labels

[indicesX, indicesY] = pol2cart([1:116]/116*2*pi, ones); %Realized you can just do plot() with 'layout','circle' argument, but this works too.

addpath('C:\Users\gsharpminor\fieldtrip-20190404');
ft_defaults
atlas = ft_read_atlas('C:\Users\gsharpminor\fieldtrip-20190404\template\atlas\aal\ROI_MNI_V4.nii'); %Tissuelabels are in atlas.tissuelabel

%Would be helpful to color the connections in the circle either by
%frontal-occipital orientation or by involvement in default mode network.
%Also might want to label the most significant parcels
%{
figure, scatter(indicesX, indicesY, 'rx'), axis equal
h = text(indicesX*1.05, indicesY*1.05, atlas.tissuelabel);
set(h, {'Rotation'}, num2cell(linspace(0, 360, 116)')); 
set(h, {'Interpreter'}, {'none'});
axis off
%}
%% 3. Load (trials)x(parcelpow) data, and calculate TLCCs between parcels

fprintf(['Loading and calculating parcelcorrs...', '\n']);
load(['D:\CU_Anschutz_SEDtoQUI\MEG_stuff\AllSubjectsAALparcelTrialPow_', band, '.mat']) %produced by: OS_MEGconvertMNEvoxelpowersToAtlas.m

parcelcorrs = [];

for k = 1:length(trialParcelPows)
    %fprintf([num2str(k), '\n']);
    tlcc = corr(trialParcelPows{k}(1:end-dt, :), trialParcelPows{k}(1+dt:end, :)); %        Time Lagged Cross Correlation, TLCC
        %ROW indices are unshifted timeseries ('causes'); COLUMN indices are future-shifted timeseries ('effects')
    parcelcorrs = cat(3, parcelcorrs, tlcc);
end %               could also use crosscorr()

%% 4. MST between parcels, all patients

fprintf(['Calculating MSTs for all patients...', '\n']);
totalConnWts = [];
nodeDegs = [];
k = {'b', 'k', 'r'};
parcelcorrsSymm = zeros(size(parcelcorrs)); %                   store the symmetrized correlations (possibly with adjacency effect subtracted)
if (dt == 0 && subtractParcAdj == true) 
    load(['D:\CU_Anschutz_SEDtoQUI\MEG_stuff\AALatlas_ParcelAdjacencies.mat']) %load AAL adjacency matrix
end

figure
for u = 1:length(trialParcelPows)
    C = abs(parcelcorrs(:, :, u)); %                    Take absolute value of TLCCs BEFORE getting max or averaging
    causalitySign = sign(C - C');  %                   '1' indicates i --> j causality; '-1' indicates j --> i causality
    if strcmp(symmetrizeMethod, 'max')
        C2 = (C > C').*C; %                                          between ij and ji, take the greater; zero out the other
        C = C2 + C2' + (C == C').*C; %                         Where ij > ji, copy ij to ji; where ij = ji, just leave them as they were.
    elseif strcmp(symmetrizeMethod, 'avg') 
        C = (C + C')/2; %        Off-diagonal mirror-pairs are averaged; diagonals just stay the same
    end
    
    if (dt == 0 && subtractParcAdj == true) %                                    subtract out the parcel adjacency effects for dt=0
        corr_adjDiff = sum(sum(C.*adjacency))/sum(sum(adjacency))...
            - sum(sum(C.*(1-adjacency)))/sum(sum(1-adjacency)); %           Get the average tlcc difference between adjacent and non-adjacent parcels 
        C = C - adjacency*(corr_adjDiff); %                                             Subtract this difference from the adjacent correlations in C
    end
    parcelcorrsSymm(:, :, u) = C;
    
    G = graph(1-C, 'upper'); % Want *greater* absolute correlation to mean lesser weight/distance & hence greater chance of being in the MST, so subtract abs(correlation) from 1.
    parcelMST = minspantree(G);
    
    subplot(8, 12, u)
    mplot = plot(parcelMST,'XData', indicesX, 'YData',indicesY, 'EdgeColor',k{labelvecNMP(u)}); axis equal, axis off  %Why is there k{}? k is not a cell object!
    h = title(names{u});
    h.Interpreter = 'none';

    %MST stats
    totalConnWts = [totalConnWts, parcelMST.Edges.Weight]; %                list of all connection weights for all subject MSTs (sec 5, below)
    nodeDegs = [nodeDegs, degree(parcelMST)]; %                                 list of all node degrees for all subject MSTs (sec 5, below)
end

suptitle({['MSTs, all patients, ', band, ', time offset = ', num2str(dt)], ['blue = Normal, black = MCI, red = PDD']})

%Dice similarity? Intersection over union? Probability of n edges being the
%same based on Cayley's Theorem? Percent difference?

%% 5. MST properties aggregated by class: edge lengths and node degrees

disp('Edge Lengths KS test------------------------------------')
%figure
for k = 1:3
%     subplot(3, 1, k)
    u = totalConnWts(:, indicesNMP{k});
    v = totalConnWts(:, indicesNMP{mod(k, 3)+1});
    fprintf([classes{k}, ' versus ', classes{mod(k, 3)+1}])
    [~, p] = kstest2(u(:), v(:)); %KS test gives incredible p-values...
    disp(p)
%     histogram(u(:), [min(totalConnWts(:)):0.01:max(totalConnWts(:))]) 
%     title([classes{k}])
end
%suptitle(['MST individual connection distances, ', band, ', time offset = ', num2str(dt)])

%Sanity check: randomly permute all distances and redo KS test
q = randperm(length(totalConnWts(:)));
u = totalConnWts(:);
u = u(q);
disp('Edge Lengths KS test, randomly permuted-------------------------------')
[~, p] = kstest2(u(1:1000), u(1001:2000));
disp(p)

disp('Degree Distribution KS test--------------------------------------')
%figure
for k = 1:3
%     subplot(3, 1, k)
    u = nodeDegs(:, indicesNMP{k});
    v = nodeDegs(:, indicesNMP{mod(k, 3)+1});
    fprintf([classes{k}, ' versus ', classes{mod(k, 3)+1}])
    [~, p] = kstest2(u(:), v(:));
    disp(p)
%     histogram(u(:), [1:1:6])
%     title([classes{k}])
end
%suptitle(['MST node degrees, ', band, ', time offset = ', num2str(dt)])

leafFrxns = sum(nodeDegs == 1, 1)/length(nodeDegs);
% figure
% bar([mean(leafFrxns(indicesNMP{1})), mean(leafFrxns(indicesNMP{2})), mean(leafFrxns(indicesNMP{3}))])
% set(gca, 'xticklabels', classes)
% title(['MST leaf fraction, ', band, ',  time offset = ', num2str(dt)])

%One-way ANOVA to see if leaf fraction differences are significant
%anova1(leafFrxns, labelvecNMP); %red lines in anova plot are median, not mean


%% 6. Create average MST between parcels, by disease class

fprintf(['Creating average per-class MSTs...', '\n']);
parcelcorrsNMP = {[], [], []}; %        this will contain the averaged correlations for all patients in each class
clims = [];
parcelcausalNMP = {[], [], []}; %       this will contain the causality direction

for k = 1:3
    parcelcorrsNMP{k} = mean(parcelcorrsSymm(:, :, indicesNMP{k}), 3); %For averaging correlations, use arithmetic mean? Harmonic? Geometric? Which is "correct"? (harmonic gets used in "average shortest path length")
%     formatSpec = 'max: %6.4f    min: %5.2f     above 0.3: %5.0f \n';
%     fprintf(formatSpec, max(parcelcorrsNMP{k}(:)), min(parcelcorrsNMP{k}(:)), sum(sum(parcelcorrsNMP{k} >= 0.3)))
    clims = [clims; [min(parcelcorrsNMP{k}(:)), max(parcelcorrsNMP{k}(:))]];
    parcelcausalNMP{k} = parcelcorrsNMP{k} - parcelcorrsNMP{k}'; % Matrix of ij - ji values. Should be positive if the row parcel mainly affects the column parcel, negative if vice versa
end

load(['D:\CU_Anschutz_SEDtoQUI\MEG_stuff\AALatlas_ParcelCentroids.mat']) %load parcel centroids for 3D figure
figure%(length(findobj('Type', 'Figure'))+1)

for u = 1:3
    G = graph(1-parcelcorrsNMP{u}, 'upper'); 
    parcelMST = minspantree(G);
    
    subplot(3, 3, 3*u-2);
    mplot = plot(parcelMST, 'layout','force'); axis equal 
    %labelnode(mplot, 1:116, atlas.tissuelabel); 
    %title([classes{u}]);
    jet(116);
    for k = 1:116
        highlight(mplot, [k], 'NodeColor',ans(k, :));
    end

    subplot(3, 3, 3*u-1);
    plot(parcelMST, 'XData',indicesX, 'YData',indicesY); axis equal, axis off  
    title([classes{u}]);
    
    subplot(3, 3, 3*u); %For dt~=0, this should really be a digraph(), with arrows showing causality
    LWidths = (1-parcelMST.Edges.Weight)/min(1-parcelMST.Edges.Weight); 
    LWidths = (LWidths - 0.95)*5;
    mplot = plot(parcelMST, 'Xdata',centroids(:, 1), 'Ydata',centroids(:, 2), 'Zdata',centroids(:, 3), 'LineWidth',LWidths, 'MarkerSize',4); axis vis3d, colormap jet
    jet(116);
    for k = 1:116
        highlight(mplot, [k], 'NodeColor',ans(k, :));
    end

end

suptitle(['MST, averaged correlation by class, ', band, ', time offset = ', num2str(dt)]);

%% 6.1. Display class-average pairwise correlation maps

figure
for k = 1:3
    subplot(1, 3, k), imagesc(parcelcorrsNMP{k}, [min(clims(:)), max(clims(:))]), axis equal, axis tight
    title(classes{k})
    %subplot(2, 3, k+3), histogram(abs(1-parcelcorrsNMP{k}))
end

cbax = axes('visible', 'off', 'Position',[0.05 0.1 .96 0.8]);
caxis(cbax, [min(clims(:)), max(clims(:))]);
h = colorbar(cbax, 'location', 'east', 'position', [0 0 1 1], 'AxisLocation','in');
suptitle(['Inter-parcel correlation map, ', band, ', time offset = ', num2str(dt)])

%% 6.2. Display class-average pairwise correlation variability (stdev) across patients 

clims = [];
h = {};
for k = 1:3
    h{k} = std(parcelcorrs(:, :, indicesNMP{k}), 0, 3);
    clims = [clims; [min(h{k}(:)), max(h{k}(:))]];
end
clims = [min(clims(:)), max(clims(:))];

figure
for k = 1:3
    subplot(1, 3, k), imagesc(h{k}, clims), axis image
    title(classes{k})
end

cbax = axes('visible', 'off', 'Position',[0.05 0.1 .96 0.8]);
caxis(cbax, clims);
colorbar(cbax, 'location', 'east', 'position', [0 0 1 1], 'AxisLocation','in');

suptitle(['Correlation variability between parcels (std. dev.), ', band, ', time offset = ', num2str(dt)])

%% 6.3. Display the class average circular MSTs with parcel labels

% This is basically a repetition of the MST display in 6, but with parcel
% labels
figure
for k = 1:3
    G = graph(1-abs(parcelcorrsNMP{k}), 'upper');
    parcelMST = minspantree(G);
    
    subplot(1, 3, k), plot(parcelMST, 'XData',indicesX, 'YData',indicesY); axis equal, axis off
    h = text(indicesX*1.06, indicesY*1.06, atlas.tissuelabel, 'FontSize',5);
    set(h, {'Rotation'}, num2cell(linspace(0, 360, 116)'));
    set(h, {'Interpreter'}, {'none'});
    axis off
    title([classes{k}]);
end

suptitle(['MST, averaged correlation by class, ', band, ', time offset = ', num2str(dt)]);

%% 7. Create class-average k-MSTs 

fprintf(['Creating class-average kMSTs...', '\n']);
kMSTdegrees = {};
kMSTedgesAll = {};
figure
for u = 1:3
    G = graph(1-abs(parcelcorrsNMP{u}), 'upper');  %create full graph
    MSTedges = [];
    for g = 1:5
        parcelMST = minspantree(G);
        q = ismember(G.Edges, parcelMST.Edges); %Find any edges already in the previous MSTs
        r = 1:length(G.Edges{:, 2});
        G = rmedge(G, r(q)); %Remove the MST edges from the full graph
        MSTedges = [MSTedges; parcelMST.Edges]; %Add the edges from the g-th MST
    end
        % Now plot the combined MSTs... specifying edge lengths changes nothing
    G = graph(MSTedges{:, 1}(:, 1), MSTedges{:, 1}(:, 2));
    subplot(1, 3, u), mplot = plot(G, 'layout','force'); axis image %The real trick is in the way you transform the k-MST onto 2D representation.
    title([classes{u}]);
    kMSTedgesAll{u} = MSTedges;
	kMSTdegrees{u} = degree(G); 
end
suptitle('kMST for the three cog classes')

for k = 1:3
    [dg, idx] = sort(kMSTdegrees{k}, 'descend');
    parcel = atlas.tissuelabel(idx(1:8))';
    degrees = dg(1:8);
%     fprintf([classes{k}, '\n'])
%     table(parcel, degrees)
%     figure, histogram(kMSTdegrees{k}, [1:1:120])
end

%% 8. Calculate shortest path length, clustering coefficient

fprintf(['Calculating shortest path length and clustering coefficients...', '\n']);
% degree distribution (already did this in section 5)
% figure
% for r = 1:3
%     histogram(kMSTdegrees{r}, [1:5:120+0])
%     hold on;
% end

% average shortest path length, L 
%   (harmonic mean of shortest distances between all node pairs in the MST. Use shortestpath(). )
%   (this can be done in a 1-MST, and actually would be simpler since in a 1-MST there is only one path between 2 nodes); 
shortestPathsNMP = cell(1, 3);
meanShortestdist = [];
xlim = 0;
for r = 1:3
    G = graph(1-abs(parcelcorrsNMP{r}), 'upper'); %         don't really need the abs(), it was done above
    parcelMST = minspantree(G);
    shortestdists = zeros(116, 116);
    for k = 1:116
        for u = 1:116
            [~, shortestdists(k, u)] = shortestpath(parcelMST, k, u); 
        end
    end 
    meanShortestdist(r) = mean(shortestdists(:));
    shortestPathsNMP{r} = shortestdists;
    if max(shortestdists(:)) > xlim %           for figure limits
        xlim = max(shortestdists(:));
    end
end
figure
fprintf(['>>>>shortest paths, KS test, NM MP PN...', '\n']);
for r = 1:3
    k = {[0 0.3 0.8], [1 0.5 0], [1 0.8 0]};
    subplot(3, 1, r), histogram(shortestPathsNMP{r}(:), [0:1:ceil(xlim)], 'FaceColor',k{r}); %ylim([0 3000]);
%     histogram(shortestdists(:), [r*0.3:1:15+r*0.3], 'BarWidth',0.5)
%     hold on;
    [~, p] = kstest2(shortestPathsNMP{r}(:), shortestPathsNMP{mod(r, 3)+1}(:));     % K-S testing on shortest paths 
    disp(p)
end
suptitle({['Histogram of 1-MST shortest path lengths, ', band, ', time offset = ', num2str(dt)]; ['blue = non-dementia, orange = MCI, gold = PDD']});
fprintf(['>>>>>mean shortest distances, NMP...', '\n']);
disp(meanShortestdist)

% clustering coefficient, C
%   (sum over all nodes of 2*(number of triangles at node)/(degree of node)/(degree of node - 1) )
%   (the number of triangles is just the number of pairs of neighbor nodes that are themselves connected)
kMSTclustCoeffs = zeros(116, 3);
for u = 1:3
    for k = 1:116
        %Find all the neighbors of node k
        q = [kMSTedgesAll{u}.EndNodes(kMSTedgesAll{u}.EndNodes(:, 1)==k, 2); ...
               kMSTedgesAll{u}.EndNodes(kMSTedgesAll{u}.EndNodes(:, 2)==k, 1)];
        %Find which of these neighbors are mutually connected
        %       number of triangles at the node: use ismember() to find all elements of q
        %       in EndNodes, then count only cases where both end nodes are in q
        %       (hence row sum is 2). 
        nTriang = sum(sum(ismember(kMSTedgesAll{u}.EndNodes, q), 2) == 2); 
        dg = length(q); %degree is just number of neighbors, hence elements in q
        kMSTclustCoeffs(k, u) = 2*(nTriang)/dg/(dg - 1); %      Basically, the fraction of edge pairs at node k that are connected into a triangle
    end
end
figure
fprintf(['>>>>clust. coeffs., KS test, NM MP PN...', '\n']);
for u = 1:3
    k = {[0 0.3 0.8], [1 0.5 0], [1 0.8 0]};
    subplot(3, 1, u), histogram(kMSTclustCoeffs(:, u), [0:0.06:max(kMSTclustCoeffs(:))], 'FaceColor',k{u}); %ylim([0 30]);
%     histogram(kMSTclustCoeffs(:, u), [u*0.01-0.06:0.06:1+u*0.01])
%     hold on;
    %NOTE: put in K-S testing between the 3 classes to quantify significance
    [~, p] = kstest2(kMSTclustCoeffs(:, u), kMSTclustCoeffs(:, mod(u, 3)+1));     % K-S testing on clust. coeffs
    disp(p)
end
suptitle({['Histogram of 5-MST node clustering coefficients, ', band, ', time offset = ', num2str(dt)]; ['blue = non-dementia, orange = MCI, gold = PDD']});
fprintf(['>>>>>mean clustering coefficients, NMP...', '\n']);
disp(sum(kMSTclustCoeffs, 1))
    
% small-worldness, S
%   (ratio of C of the network to C of a random graph, divided by the ratio of L of the network to L of a random graph)
%   (conclude small-worldness when S > 1)

% modularity, Q
%   (calculating this needs pre-defined 'modules' of nodes. Don't do it)

% then significance testing.
% Could do this for individual patients, but class average kMST is simpler for now

data = struct(); %Save this stuff as a structure object
data.shortestPaths_NMP = shortestPathsNMP;
data.kMSTclustCoeffs_NMP = kMSTclustCoeffs;

%% 9. Create a distance matrix calculated by using the correlations as coordinates, then generate MST, cut at hybrid edges & do multivariate runs test

%Use the class-averaged correlation matrix in other words, not the
%concatenation of all individual matrices
fprintf(['Creating pooled distance matrix and 2-class MVR cut...', '\n']);

corrVects = []; %Use rows, meaning the 'present' parcel correlated to the 'future' parcel
for k = 1:3
    corrVects = [corrVects; parcelcorrsNMP{k}]; %matrix containing all the correlation vectors for all parcels for all patients (89*116=10324)
end

distances = zeros(length(corrVects)); 
for d = 1:116
    distances = distances + (corrVects(:, d) - corrVects(:, d)').^2; 
end
distances = sqrt(distances);
%figure, imagesc(distances)
%distances = randn(length(corrVects)); %sanity check: random normally distributed points should yield result close to expectation

for k = 1:3
    u = [k, mod(k, 3)+1];
    v = [u(1)*116-115:u(1)*116, u(2)*116-115:u(2)*116];
    %v = [indicesNMP{u(1)}(1)*116-115:indicesNMP{u(1)}(end)*116, indicesNMP{u(2)}(1)*116-115:indicesNMP{u(2)}(end)*116];
    G = graph(distances(v, v), 'upper');  %Take only the parcel distance pairs between the classes being compared
    parcelMST = minspantree(G);
    C = sum(degree(parcelMST).*(degree(parcelMST) - 1)/2); 
        % get "Function 'subsindex' is not defined for values of class 'graph'" because a variable named 'degree' somehow already exists
    
    figure 
    subplot(1, 2, 1) 
    mplot = plot(parcelMST, 'layout','force'); axis image
    highlight(mplot, [1:116],'NodeColor','b');
    highlight(mplot, [116+1:2*116],'NodeColor','r')
    
    q = ismember(parcelMST.Edges.EndNodes, 1:116); %Identify the edges that are from class k
    w = parcelMST.Edges.EndNodes(q(:, 1) ~= q(:, 2), :); %Find the 'hybrid' edges (dissimilar node classes)
    parcelMST = rmedge(parcelMST, w(:, 1), w(:, 2)); %Remove the hybrid edges from the MST

    subplot(1, 2, 2) 
    mplot = plot(parcelMST, 'layout','force'); axis image
    highlight(mplot, [1:116],'NodeColor','b');
    highlight(mplot, [116+1:2*116],'NodeColor','r')
    
    suptitle({['Class average correlation-distance MST, cut by dissimilarity, ', band, ', time offset = ', num2str(dt), ':'], ['(', classes{k}, ', blue; ', classes{mod(k, 3)+1}, ', red)']})
        
    % Then apply multidimensional two-sample test from Friedman & Rafsky 1979 (or radial Smirnov) on the resulting
    % fragmented structure. Basically, the more the MST is fragmented
    % by removing 'hybrid' connections, the more likely that the two MSTs are from
    % the same distribution. 

    m = 116;
    n = length(v) - m;
    N = n + m;

    runstestMean = 2*m*n/N + 1;
    runstestVar = 2*m*n/(N*(N - 1))*(...
        (2*m*n - N)/N + ...
            (C - N + 2)/((N - 2)*(N - 3))*...
                (N*(N - 1) - 4*m*n + 2)...
         ); 
    sqrt(runstestVar);
    
    fprintf(['>>>>Disease classes compared: ', classes{k}, ' and ', classes{mod(k, 3)+1}, '\n'])
    fprintf(['        Number of trees from hybrid edge cut (actual): ', num2str(length(w)+1), '\n'])
    fprintf(['        Number of trees from hybrid edge cut (expect): ', num2str(runstestMean), '\n'])
    fprintf(['        Standard devs from expectation: ', num2str((runstestMean - length(w) - 1)/sqrt(runstestVar)), '\n'])
    
    data.numTreesNM_MP_NP{1, k} = length(w)+1; %save the tree count
    data.numTreesNM_MP_NP{2, k} = runstestMean; %save expected # trees
    data.numTreesNM_MP_NP{3, k} = ( length(w) + 1 - runstestMean )/sqrt(runstestVar); %save the stdevs from expected
end
    
%% 10. Multivariate runs for all 3 classes at once

fprintf(['Creating pooled distance matrix and 3-class MVR cut...', '\n']);

corrVects = []; %Use rows, meaning the 'present' parcel correlated to the 'future' parcel
for k = 1:3
    corrVects = [corrVects; parcelcorrsNMP{k}]; %matrix containing all the correlation vectors for all parcels for all patients (89*116=10324)
end

distances = zeros(length(corrVects)); 
for d = 1:116
    distances = distances + (corrVects(:, d) - corrVects(:, d)').^2; 
end
distances = sqrt(distances);
%distances = randn(length(corrVects)); %sanity check: random normally distributed points should yield result close to expectation

v = [1:3*116];
G = graph(distances(v, v));  %Take only the parcel distance pairs between the classes being compared
parcelMST = minspantree(G);
MST = parcelMST; %              save a copy of the pooled MST without edges cut, for variance simulation later

figure 
subplot(1, 2, 1) 
mplot = plot(parcelMST, 'layout','force'); axis image
highlight(mplot, [1:116],'NodeColor','r');
highlight(mplot, [116+1:2*116],'NodeColor','g')
highlight(mplot, [2*116+1:3*116],'NodeColor','b')

qN = ismember(parcelMST.Edges.EndNodes, 1:116); %Identify the nodes that are from class N
qM = ismember(parcelMST.Edges.EndNodes, 116+1:2*116)*2; %Identify the nodes that are from class M
qP =  ismember(parcelMST.Edges.EndNodes, 2*116+1:3*116)*3; %Identify the nodes that are from class P

q = qN + qM + qP; %             put all the node identities together

classAdjacency = zeros(3); %create class adjacency matrix from q
for h = 1:length(q)
    classAdjacency(q(h, 1), q(h, 2)) = classAdjacency(q(h, 1), q(h, 2)) + 1;
end

classAdjacency = classAdjacency + classAdjacency' - diag(diag(classAdjacency)); %that should do it

hybridEdges = (q(:, 1) ~= q(:, 2));
w = parcelMST.Edges.EndNodes(hybridEdges, :); %Find the 'hybrid' edges (dissimilar node classes)
parcelMST = rmedge(parcelMST, w(:, 1), w(:, 2)); %Remove the hybrid edges from the MST

subplot(1, 2, 2) 
mplot = plot(parcelMST, 'layout','force'); axis image
highlight(mplot, [1:116],'NodeColor','r');
highlight(mplot, [116+1:2*116],'NodeColor','g')
highlight(mplot, [2*116+1:3*116],'NodeColor','b')

suptitle({['Class average correlation-distance MST, 3-way cut, ', band, ', time offset = ', num2str(dt), ':'], ['(', classes{1}, ', red; ', classes{2}, ', green; ', classes{3}, ', blue)']})

fprintf(['        Number of trees from hybrid edge cut (actual): ', num2str(length(w)+1), '\n'])

%...but how many of each tree color? 
treeCountsNMP = [];
conn = conncomp(parcelMST);
treeCountsNMP(1) = length(unique(conn(1:116)));
treeCountsNMP(2) = length(unique(conn(116+1:2*116)));
treeCountsNMP(3) = length(unique(conn(2*116+1:3*116)));

%data.numTrees3Classes = length(w)+1; %save the tree count
data.classAdjacencyNMP = classAdjacency; %save the class adjacency
data.treeCountsNMP = treeCountsNMP; %save the tree counts (sum gives total tree count)

%% 11. Run the kMST force-directed dimension reduction
%{
[positions, h] ...
        = OS_kMSTmaker_forceDirectedDimReduce_noRepel(corrVects, 1600, 0.015, 0.042); %Results are amazingly consistent, given it's randomly initialized!
h.CData = [repmat([0 0.3 0.8], 116, 1); repmat([1 0.8 0.1], 116, 1); repmat([1 0.2 0], 116, 1)]; %            Hmm, doesn't seem to find much separation between the 3 classes, though.
%}
%% 12. Run the node class permutation on the 3-class pooled MST to get class-adjacency statistics

fprintf(['Running 3-class node permutation on pooled MST...', '\n']);

gcp;
[data.classAdj_means, data.classAdj_stdevs, totalTrees3class_mean, totalTrees3class_stdev]...
        = OS_MultivariateRuns_varianceByPermuting(MST, [116 116 116], 50000);

fprintf(['>>>>Class adjacency matrix:', '\n']);
disp(data.classAdjacencyNMP)
fprintf(['>>>>Class adjacency deviations from expected:', '\n']);
disp((data.classAdjacencyNMP - data.classAdj_means)./data.classAdj_stdevs)
fprintf(['>>>>Total trees, 3-class: expected, stdev, deviation from expect...', '\n']);
disp([totalTrees3class_mean, totalTrees3class_stdev, (sum(data.treeCountsNMP) - totalTrees3class_mean)/totalTrees3class_stdev])

%% 13. Save the data structure

destpath = strcat('D:\CU_Anschutz_SEDtoQUI\MEG_stuff\MEG_paper!\Paper2_connectivity\PAPER2_connectivity_figures\201208MSTresults\201208MSTresults_', band, '_dt', num2str(dt), '.mat');

save(destpath, 'data')

toc


