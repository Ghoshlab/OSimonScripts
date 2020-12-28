% 1. Load all the sourceratio voxels for all patients
% 2. Remove all NaN entries from each patient
% 3. Run MANOVA on the three disease groups 

%clear;
addpath('C:\Users\gsharpminor\fieldtrip-20190404');
ft_defaults
bandnames = {'alpha', 'beta', 'delta', 'gamma', 'theta'};
[b, ~] = listdlg('PromptString','Select brainwave band:', 'SelectionMode','single', 'ListString',bandnames);
band = bandnames{b};
destdir = 'G:\201124shuttlevector\';
template_mri = ft_read_mri('C:\Users\gsharpminor\fieldtrip-20190404\external\spm8\templates\T1.nii');


%% Get the parcellation set up

atlasChoice = 1;

%Load the source ratio data
load('D:\MEGdata_ParkinsonsMCI\BK_K02_Normal_cases\BK_K02_11\BK_K02_11_sources.mat','intalpharatio'); 

%Load the nice & detailed 2016 Brainnetome Atlas...
%atlas = ft_read_atlas('D:\MEGdata_ParkinsonsMCI\templates\BN_Atlas_246_2mm.nii'); %OS-- Brainnetome atlas
%atlas = ft_read_atlas('C:\Users\gsharpminor\fieldtrip-20190404\template\atlas\brainnetome\BNA_MPM_thr25_1.25mm.nii'); OS-- this one gives 'index exceeds matrix dimensions'. Unsure why

%Or there is also the less detailed AAL atlas... 
atlases = {'C:\Users\gsharpminor\fieldtrip-20190404\template\atlas\aal\ROI_MNI_V4.nii',...
    'C:\Users\gsharpminor\fieldtrip-20190404\template\atlas\brodmann.nii',...
    'C:\Users\gsharpminor\fieldtrip-20190404\template\atlas\labels_10July.nii'};
atlas = ft_read_atlas(atlases{atlasChoice});

if atlasChoice == 1

    % parcellate the sources into atlas subdivisions
    cfg=[];
    cfg.method = 'eig';
    cfg.parcellation = 'tissue';
    cfg.parameter='pow';
    atlas.pos = intalpharatio.pos; % otherwise error not in same space
    intalpharatioparc = ft_sourceparcellate(cfg,intalpharatio,atlas); 

    % plot parcellation
    cfg=[];
    cfg.parameter='pow';
    cfg.method='ortho';
    cfg.funparameter='pow';
    %cfg.atlas = atlas;
    ft_sourceplot(cfg,intalpharatioparc,atlas);
    colormap jet;

    %clear('intalpharatio')

 elseif atlasChoice == 2
     
    %The size of the interpolated data is 91x109x91, while this parcellation is 181x217x181. 
    %So the parcellation has to be warped/resized to fit
    t =...
        [1/2          0             0             0;
         0             1/2          0             0;
         0             0             1/2          0;
         0             0             0             1];
     
    tform = affine3d(t);
    atlas.brick0 = imwarp(atlas.brick0, tform); 
    atlas.dim = [91 109 91];
    
    %parcellate the sources into atlas subdivisions
    cfg=[];
    cfg.method = 'eig';
    cfg.parcellation = 'brick0';
    cfg.parameter='pow';
    atlas.pos = intalpharatio.pos; % otherwise error not in same space
    intalpharatioparc = ft_sourceparcellate(cfg, intalpharatio, atlas); 

    % plot parcellation
    cfg=[];
    cfg.parameter='pow';
    cfg.method='ortho';
    cfg.funparameter='brick0';
    %cfg.atlas = atlas; 
    ft_sourceplot(cfg, intalpharatioparc, atlas), colormap jet;
    
    %rename the fields to go with fieldnames from AAL
    atlas.tissue = atlas.brick0;
    atlas.tissuelabel = atlas.brick0label;
    atlas = rmfield(atlas, {'brick0', 'brick0label'});
    
elseif atlasChoice == 3
     
    %The size of the interpolated data is 91x109x91, while this parcellation is 182x218x182. 
    %So the parcellation has to be warped/resized to fit

    t =...
        [1/2          0             0             0;
         0             1/2          0             0;
         0             0             1/2          0;
         0             0             0             1];
     
    tform = affine3d(t);
    atlas.brick0 = imwarp(atlas.brick0, tform); 
    atlas.dim = [91 109 91];
    
    %parcellate the sources into atlas subdivisions
    cfg=[];
    cfg.method = 'eig';
    cfg.parcellation = 'brick0';
    cfg.parameter='pow';
    atlas.pos = intalpharatio.pos; % otherwise error not in same space
    intalpharatioparc = ft_sourceparcellate(cfg, intalpharatio, atlas); 

    % plot parcellation
    cfg=[];
    cfg.parameter='pow';
    cfg.method='ortho';
    cfg.funparameter='brick0';
    %cfg.atlas = atlas; 
    ft_sourceplot(cfg, intalpharatioparc, atlas), colormap jet;
    
    %rename the fields to go with fieldnames from AAL
    atlas.tissue = atlas.brick0;
    atlas.tissuelabel = atlas.brick0label;
    atlas = rmfield(atlas, {'brick0', 'brick0label'});
    
end

atlascopy = zeros(size(atlas.tissue)); % this will hold the neg-log probs


%% Get all the sourceratio voxels for all patients (basically same code as in OS_MEGrandomForest01)

c1 = {'BK_K02_Normal', 'BK_KL2_Normal', 'BK_K02_MCI', 'R21_PDD'};
sourceratios = {};
X = []; %observation vector (dependent variables)
Y = []; %labels vector (independent variable)

for r = 1:length(c1) %Loop through all 3 conditions
    a = dir(['D:\MEGdata_ParkinsonsMCI\', c1{r}, '_cases\']);
    for s = 3:length(a) %Get all the sources.mat files for the disease class
        a(s).name
        load([a(s).folder, '\', a(s).name, '\', a(s).name, '_sources.mat'], ['int', band, 'ratio']); %Don't think the interpolated source space would add any useful info
        sourceratio = eval(['int', band, 'ratio']);
        A = zeros(sourceratio.dim);
        A(:) = sourceratio.pow;
        sourceratios = [sourceratios, A]; %no real need to keep them as 3D arrays
        Y = [Y, r]; %class labels
    end
end

Y(Y == 1) = 2; %So all Normal patients have same label

eval(['sourceratios_', band, '= sourceratios']); %keep it for later


%% Subset by each parcel, turning the ratio voxels into average intensities per parcel

nParcels = length(atlas.tissuelabel);
nSubjects = length(Y);
subjParcel = zeros(nSubjects, nParcels);

for x = 1:nParcels
    parcel = (atlas.tissue == x); %Make sure the atlas and the source data are not flipped relative to each other... (done)
    runsum = [];
    for y = 1:nSubjects
        runsum = [runsum, sum(sourceratios{y}(parcel))]; %for each subject sum the ratios within the parcel
    end 
    subjParcel(:, x) = runsum/sum(parcel(:)); %Divide each entry by the number of voxels in the corresponding parcel
end


%% Carry out ANOVA for each parcel across subjects

for q = 1:nParcels
   [p, ~, stats] = anova1(subjParcel(:, q), Y, 'off'); %need to store the multiple comparisons data too
    parcels_p(q) = p;
    parcels_stats{q} = stats;
end
indices = find(parcels_p<0.01); 

%list the names of the parcels where p<0.01
atlas.tissuelabel(indices)
parcels_p(indices)

% % Sanity-check: randomly permute class labels and see if significance drops (it should)
% permute = randperm(numel(Y)); Yperm = Y(permute); %randomly shuffle the positions
% 
% for q = 1:nParcels
%     parcels_p_perm(q) = anova1(subjParcel(:, q), Yperm, 'off'); %need to store the multiple comparisons data too
% end
% 
% atlas.tissuelabel(parcels_p_perm<0.01)

p_output = table(atlas.tissuelabel', parcels_p', [1:116]', 'VariableNames',{'parcel_name' 'parcel_p_value' 'parcel_ID'});
%writetable(p_output, [destdir, 'parcelPvalues_', band, '.csv'])

%p_output = table(atlas.tissuelabel(indices)', parcels_p(indices)');
%writetable(p_output, [destdir, 'parcelPvalues0.01_', band, '.csv'])

[x, q] = sort(parcels_p');
% parcel_name = atlas.tissuelabel(q)'; p_val = parcels_p(q)'; %Huh? Why make the same table a third time? Oh, I was saving them...
% p_output = table(parcel_name, p_val);
%writetable(p_output, [destdir, 'parcelPvalues_', band, '_sorted.csv'])

%% Create a line plot showing the means for NMP for the six most significant parcels

n = [];

for u = 1:6
    n = [n; parcels_stats{q(u)}.means];
end
p = {'Normal', 'MCI', 'PDD'};
figure, plot(n')
set(gca, 'xtick', [1:3], 'xticklabel',p)
title([band, ' ratio power as function of disease class, top 6 parcels'])
%legend(m(1:6), 'Interpreter','none') %PROBLEM: why is it now stopping at u = 2 ???
%savefig([destdir, 'lineplotMostSigParcels_', band, '.fig'])


%% View the multiple comparisons plot for each significant parcel

parcels_stats_sig = parcels_stats(indices); %Just want the stats from the significant cases
figure
for n = 1:4%length(parcels_stats)
    %subplot(6, 6, n) %It won't let me put them in a subplot!
    multcompare(parcels_stats_sig{n}) %Why does multcompare() not do anything if I run it a second time?
        %Could just save the means and error bars and plot
        %separately...
    title(atlas.tissuelabel(indices(n)), 'Interpreter','none')
    h = msgbox('Next one?');
    uiwait(h);
end %It is showing *decreasing power for the lower frequencies in MCI & PDD, not increasing, because this is RATIO power


 %% ANOVA of whole-brain power

sourceratios_mean = [];
for u = 1:length(sourceratios)
    sourceratios_mean = [sourceratios_mean, mean(sourceratios{u}(sourceratio.inside))];
end

anova1(-sourceratios_mean, Y) %Where is the change of sign coming from?


%% Benjamini-Hochberg procedure, for FDR control
    %arrange parcel P-values in increasing order, then take all values
    %below the largest k such that P(k)<=pval/116*k.
alpha = 0.05; %only alpha and beta parcels remain at 0.04. Almost half the parcels are significant in beta.
s = sortrows(p_output, 'parcel_p_value');
m = max(find(s{:, 2}' <= alpha/116*(1:116))); %[~, m] = max()
figure, plot(1:116, s{:, 2})
hold on;
plot(1:116, alpha/116*(1:116), 'r') 
title({['Benjamini-Hochberg result for parcel activation differences:'], [band, ', red line for p<', num2str(alpha)]})
p_bh = s(1:m, :) 
% writetable(p_bh, [destdir, 'parcelPvaluesFDR_', band, '.csv'])
% writetable(s, [destdir, 'parcelPvalues_', band, '.csv'])


%% Color the parcels by significance and display: version iii (like i, but with MRI anatomy added for context)

for x = 1:m %use Benjamani-Hochberg parcels only
    parcel = (atlas.tissue == p_bh{x, 3});
    atlascopy(parcel) = -log10(p_bh{x, 2}); 
end

clims = [min(-log10(p_bh{:, 2})), max(-log10(p_bh{:, 2}))];
if isempty(clims)
    clims = [1, 2];
elseif clims(1) == clims(2)
    clims = [clims(1)-1, clims(2)+1];
end

figure
for u = 1:25
    slice = squeeze(template_mri.anatomy(85-3*u, :, :));
    %make a mask to remove the parts of the first that overlap with the second, then sum them together.
    probs = squeeze(atlascopy(85-3*u, :, :));
    thresh = clims(1);
    sliceMasked = slice.*(probs<thresh); %mask out the part of MRI that overlaps active parcels
    %probsMasked = probs.*(probs>thresh);
    pmax = clims(2); %max(probs(:));
    map = jet(256);
    probsColorized = ind2rgb(uint8((probs-thresh)*256/(pmax-thresh)), map);
    probsColorized(:, :, 3) = probsColorized(:, :, 3).*(probsColorized(:, :, 3)>map(1,3)); %get rid of the background blue
    subplot(5, 5, u); imshow(imrotate(probsColorized + repmat(sliceMasked, 1, 1, 3), 90)) 
    colormap jet; 
end 
cbax = axes('visible', 'off', 'Position',[0.05 0.1 .96 0.8]); 
caxis(cbax, clims);
h = colorbar(cbax, 'location', 'east', 'position', [0 0 1 1], 'AxisLocation','in');
h = suptitle(['Negative log of significance: ', band, 'ratio (left to right)']);
h.Interpreter = 'none';
%saveas(gcf, [destdir, 'Figure2_FDR', num2str(alpha), '_', band, '.jpg'])

%{
%% Color the parcels by significance and display: version ii, "glass brain"

for x = 1:m %use Benjamani-Hochberg parcels only
    parcel = (atlas.tissue == p_bh{x, 3});
    atlascopy(parcel) = -log10(p_bh{x, 2}); 
end

%{
New display method: 
1) add together neg-log significances through X, Y and Z directions
2) project these significances onto a midway slice of the template MRI in
the YZ, XZ, and XY planes
%}

figure

% XY plane
midslice = template_mri.anatomy(:, :, 45);
%make a mask to remove the parts of the first that overlap with the second, then sum them together.
probsFlattened = squeeze(sum(atlascopy, 3));
thresh = 40; %thresholding here is kind of arbitrary... want to show the structure of the significant regions as well as possible...
midsliceMasked = midslice.*(probsFlattened<thresh); 
probsMasked = probsFlattened.*(probsFlattened>thresh);
pmax = max(probsFlattened(:));
map = jet(256);
probsColorized = ind2rgb(uint8((probsFlattened-thresh)*256/(pmax-thresh)), map);
probsColorized(:, :, 3) = probsColorized(:, :, 3).*(probsColorized(:, :, 3)>map(1,3)); %get rid of the background blue
subplot(1, 3, 1); imshow(imrotate(probsColorized + repmat(midsliceMasked, 1, 1, 3), 90)) %Does it need a colorbar? Somewhat arbitrary, since we are looking through-space


% XZ plane
midslice = squeeze(template_mri.anatomy(:, 54, :));
probsFlattened = squeeze(sum(atlascopy, 2));
thresh = 40; 
midsliceMasked = midslice.*(probsFlattened<thresh); 
probsMasked = probsFlattened.*(probsFlattened>thresh);
pmax = max(probsFlattened(:));
map = jet(256);
probsColorized = ind2rgb(uint8((probsFlattened-thresh)*256/(pmax-thresh)), map);
probsColorized(:, :, 3) = probsColorized(:, :, 3).*(probsColorized(:, :, 3)>map(1,3)); 
subplot(1, 3, 2); imshow(imrotate(probsColorized + repmat(midsliceMasked, 1, 1, 3), 90))


% YZ plane
midslice = squeeze(template_mri.anatomy(45, :, :));
probsFlattened = squeeze(sum(atlascopy, 1));
thresh = 40; 
midsliceMasked = midslice.*(probsFlattened<thresh); 
probsMasked = probsFlattened.*(probsFlattened>thresh);
pmax = max(probsFlattened(:));
map = jet(256);
probsColorized = ind2rgb(uint8((probsFlattened-thresh)*256/(pmax-thresh)), map);
probsColorized(:, :, 3) = probsColorized(:, :, 3).*(probsColorized(:, :, 3)>map(1,3)); 
subplot(1, 3, 3); imshow(imrotate(probsColorized + repmat(midsliceMasked, 1, 1, 3), 90)) %why is this one displayed disproportionately small? It seems to be locking everything to 91 pixels width...


%% Color the parcels by significance and display: version i

for x = 1:nParcels
    parcel = (atlas.tissue == x);
    atlascopy(parcel) = -log10(parcels_p(x)); 
end

clims = [1 3.5];
figure
for m = 1:25
subplot(5, 5, m); imagesc(imrotate(squeeze(atlascopy(85-3*m, :, :)), 90), clims), axis equal, axis off;
colormap jet;
end 
cbax = axes('visible', 'off', 'Position',[0.05 0.1 .96 0.8]);
caxis(cbax, clims);
h = colorbar(cbax, 'location', 'east', 'position', [0 0 1 1], 'AxisLocation','in');
h = suptitle(['Negative log of significance: ', band, 'ratio (left to right)']);
h.Interpreter = 'none';
%savefig([destdir, 'significanceMapParcels_', band, '.fig'])
%}
%end

% %% 3D display of the parcels that are significant? 
% 
% % Go backwards: from the significant parcels to the indices of these
% % parcels to displaying the atlas regions colored by likelihood (1/p)
% 
% %[A, B] = zeros(intalpharatio.dim);
% 
% indices = 1:prod(intalpharatio.dim);
% A = indices(atlas.tissue(:) ~= 0);
% B = indices(intalpharatio.inside(:) ~= 0);
% 
% [x, y, z] = ind2sub(intalpharatio.dim, A);
% coordsA = [x', y', z'];
% [x, y, z] = ind2sub(intalpharatio.dim, B);
% coordsB = [x', y', z'];
% 
% min(coordsA), max(coordsA)
% min(coordsB), max(coordsB)
% sum((atlas.tissue(:) ~= 0) & (intalpharatio.inside(:) ~= 0)) %The atlas is completely inside the sourcespace. But it covers barely over half of the source voxels. Possibly because it leaves out deeper brain structures?
% 
% clear('A', 'B', 'x', 'y', 'z')
% 
% figure, scatter3(coordsA(:, 1), coordsA(:, 2), coordsA(:, 3)), axis vis3d, axis equal
% figure, scatter3(coordsB(:, 1), coordsB(:, 2), coordsB(:, 3)), axis vis3d, axis equal
% figure, imagesc(squeeze(atlas.tissue(45, :, :))), axis equal
% figure, imagesc(squeeze(intalpharatio.pow(45, :, :))), axis equal, colormap jet


