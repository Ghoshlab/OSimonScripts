%Still in sensor not source space, but as long as we average over the
%channels to get overall power density it should be OK

load('D:\MEGdata_ParkinsonsMCI\OS_channelListEtc.mat');
paths = {'D:\MEGdata_ParkinsonsMCI\BK_K02_Normal_cases\', 
    'D:\MEGdata_ParkinsonsMCI\BK_KL2_Normal_cases\', 
    'D:\MEGdata_ParkinsonsMCI\BK_K02_MCI_cases\', 
    'D:\MEGdata_ParkinsonsMCI\R21_PDD_cases\'
    };
allsubjects = cell(2, 1); %This stores the power at all the frequencies, together with class labels ('1'=Normal, '2'=MCI, '3'=PDD).

%1. Get average power as a function of frequency for each subject
%2. Interpolate so that there is a consistent number of frequencies sampled from 0-50Hz; average over trials and channels
%3. Repeat for all three groups, Normal/MCI/PDD
for k=1:4
    k
    datapows = {};
    a = dir(paths{k});
    for b = 3:length(a) %Loading in the datapows
        fprintf(['loading: ', a(b).name, '\n']);
        load(strcat(paths{k}, a(b).name, '\', a(b).name, '_spec.mat'), 'datapow');
        datapows{b-2} = datapow;
    end
    for c = 1:length(datapows) %Don't really need this second separate loop...
        if strcmp(datapows{b-2}.label{end}, 'EKG')
            datapows{c}.powspctrm = datapows{c}.powspctrm(:, 1:end-1, :); %remove the EKG channel
            fprintf('removed EKG... ')
        end
        if isequal(frequencies, datapows{c}.freq) == 0
            Q = interp1(datapows{c}.freq, permute(datapows{c}.powspctrm, [3, 1, 2]), frequencies, 'spline'); %Interpolation: 'frequencies' is the desired (target) frequencies set
            Q = permute(Q, [2, 3, 1]);
        else
            Q = datapows{c}.powspctrm; %Don't run interpolation if you don't have to
        end
        Q = mean(Q, 1); %Average over all trials
        Q = mean(Q, 2); %Average over all sensors (so we don't have to worry about which specific sensors were dropped)
        allsubjects{1} = [allsubjects{1}, squeeze(Q)];
        allsubjects{2} = [allsubjects{2}, k];
        fprintf([num2str(c), 'th done: ', num2str(length(datapows{c}.freq)), 'freqs \n']);
    end
end
clear('datapow', 'datapows')

allsubjects{2}(allsubjects{2}==2)=1; %combine the class labels for 'KL2' and 'K02'
allsubjects{2}(allsubjects{2}==3)=2;
allsubjects{2}(allsubjects{2}==4)=3;

%Should the powers be normalized for each subject? Or left as-is?

%% 1. Carry out one-way MANOVA to see if there is a significant within-band change between Normal-->MCI-->PDD

%(Used same frequency cutoffs as in Don's "script_resting_meg_spectral.m")
bandlimits = [8  12, 13  30, 4   7, 0 3.5, 31  55];

%Select just the frequencies corresponding to alpha band: 
fprintf(['A:_', num2str(sum(frequencies > bandlimits(1) & frequencies < bandlimits(2))), '_frequencies'])
Ra = allsubjects{1}(find(frequencies > bandlimits(1) & frequencies < bandlimits(2)), :)';
[~, alpha_p, alpha_stats] = manova1(Ra, allsubjects{2}, 0.05); %look for 0.05 significance

% ...beta band: 
fprintf(['B:_' num2str(sum(frequencies > bandlimits(3) & frequencies < bandlimits(4))), '_frequencies'])
Rb = allsubjects{1}(find(frequencies > bandlimits(3) & frequencies < bandlimits(4)), :)';
[~, beta_p, beta_stats] = manova1(Rb, allsubjects{2}, 0.05); %look for 0.05 significance

% ...delta band: 
fprintf(['D:_' num2str(sum(frequencies > bandlimits(7) & frequencies < bandlimits(8))), '_frequencies'])
Rd = allsubjects{1}(find(frequencies > bandlimits(7) & frequencies < bandlimits(8)), :)';
[~, delta_p, delta_stats] = manova1(Rd, allsubjects{2}, 0.05); %look for 0.05 significance

% ...gamma band: 
fprintf(['G:_' num2str(sum(frequencies > bandlimits(9) & frequencies < bandlimits(10))), '_frequencies'])
Rg = allsubjects{1}(find(frequencies > bandlimits(9) & frequencies < bandlimits(10)), :)';
[~, gamma_p, gamma_stats] = manova1(Rg, allsubjects{2}, 0.05); %look for 0.05 significance

% ...theta band
fprintf(['T:_' num2str(sum(frequencies > bandlimits(5) & frequencies < bandlimits(6))), '_frequencies'])
Rt = allsubjects{1}(find(frequencies > bandlimits(5) & frequencies < bandlimits(6)), :)';
[~, theta_p, theta_stats] = manova1(Rt, allsubjects{2}, 0.05); %look for 0.05 significance

% ...and all 99 frequencies
%[allfreq_d, allfreq_p, allfreq_stats] = manova1(allsubjects{1}', allsubjects{2}, 0.05); %Just an error: 'The within-group sum of squares and cross products matrix is singular.'

alpha_p
beta_p 
delta_p
gamma_p 
theta_p


%% 2. try manovacluster() to look at the relationships

figure, 
subplot(5, 1, 1), manovacluster(alpha_stats), title('MANOVA of bandpower: alpha')
subplot(5, 1, 2), manovacluster(beta_stats), title('MANOVA of bandpower: beta')
subplot(5, 1, 3), manovacluster(delta_stats), title('MANOVA of bandpower: delta')
subplot(5, 1, 4), manovacluster(gamma_stats), title('MANOVA of bandpower: gamma')
subplot(5, 1, 5), manovacluster(theta_stats), title('MANOVA of bandpower: theta')


%% 3. Examine the pairwise differences too

% Normal vs. MCI
%[allfreq_d_NvsM, allfreq_p_NvsM, ] = manova1(allsubjects{1}(:, 1:48)', allsubjects{2}(1:48), 0.05); %within-group sum of squares and cross products matrix is singular

idx = find(allsubjects{2} == 1 | allsubjects{2} == 2);
[d, alpha_p_NvsM, ] = manova1(Ra(idx, :), allsubjects{2}(idx), 0.05); 
[d, beta_p_NvsM, ] = manova1(Rb(idx, :), allsubjects{2}(idx), 0.05); 
[d, delta_p_NvsM, ] = manova1(Rd(idx, :), allsubjects{2}(idx), 0.05); 
[d, gamma_p_NvsM, ] = manova1(Rg(idx, :), allsubjects{2}(idx), 0.05); 
[d, theta_p_NvsM, ] = manova1(Rt(idx, :), allsubjects{2}(idx), 0.05); 

% Normal vs. PDD
%[allfreq_d_NvsP, allfreq_p_NvsP, allfreq_stats_NvsP] = manova1(allsubjects{1}(:, [1:11, 49:65])', allsubjects{2}([1:11, 49:65]), 0.05); %within-group sum of squares and cross products matrix is singular

idx = find(allsubjects{2} == 1 | allsubjects{2} == 3);
[d, alpha_p_NvsP, ] = manova1(Ra(idx, :), allsubjects{2}(idx), 0.05); 
[d, beta_p_NvsP, ] = manova1(Rb(idx, :), allsubjects{2}(idx), 0.05); %within-group sum of squares and cross products matrix is singular
[d, delta_p_NvsP, ] = manova1(Rd(idx, :), allsubjects{2}(idx), 0.05); 
[d, gamma_p_NvsP, ] = manova1(Rg(idx, :), allsubjects{2}(idx), 0.05); %within-group sum of squares and cross products matrix is singular
[d, theta_p_NvsP, ] = manova1(Rt(idx, :), allsubjects{2}(idx), 0.05); 

% MCI vs. PDD
%[allfreq_d_MvsP, allfreq_p_MvsP, allfreq_stats_MvsP] = manova1(allsubjects{1}(:, 12:65)', allsubjects{2}(12:65), 0.05); %within-group sum of squares and cross products matrix is singular

idx = find(allsubjects{2} == 2 | allsubjects{2} == 3);
[d, alpha_p_MvsP, ] = manova1(Ra(idx, :), allsubjects{2}(idx), 0.05); 
[d, beta_p_MvsP, ] = manova1(Rb(idx, :), allsubjects{2}(idx), 0.05); 
[d, delta_p_MvsP, ] = manova1(Rd(idx, :), allsubjects{2}(idx), 0.05); 
[d, gamma_p_MvsP, ] = manova1(Rg(idx, :), allsubjects{2}(idx), 0.05); 
[d, theta_p_MvsP, ] = manova1(Rt(idx, :), allsubjects{2}(idx), 0.05); 

%allfreq_p_NvsM

alpha_p_NvsM
beta_p_NvsM
delta_p_NvsM
gamma_p_NvsM
theta_p_NvsM

%allfreq_p_NvsP

alpha_p_NvsP
beta_p_NvsP
delta_p_NvsP
gamma_p_NvsP
theta_p_NvsP

%allfreq_p_MvsP

alpha_p_MvsP
beta_p_MvsP
delta_p_MvsP
gamma_p_MvsP
theta_p_MvsP

%Any 'matrix is singular' errors are usually because the number of observations
%(frequencies) is larger than the total number of subjects.

%Essentially, any significant effect found here is due to the fact that the
%matrix is close to singular, unless you narrow it down to a single band.

%Try also between-band correlations, and see if those differ between
%groups.


%% 4. Comparison of total all-band power across subjects with ANOVA

RaTot = mean(Ra, 2); %Use mean not sum, as different bands don't have the same number of frequencies
RbTot = mean(Rb, 2);
RdTot = mean(Rd, 2);
RgTot = mean(Rg, 2);
RtTot = mean(Rt, 2);

alphapower_p = anova1(RaTot, allsubjects{2}) %p < 0.1205 n.s.
title('ANOVA: alpha power')
xticklabels({'Normcog', 'MCI', 'PDD'})
betapower_p = anova1(RbTot, allsubjects{2})
title('ANOVA: beta power')
xticklabels({'Normcog', 'MCI', 'PDD'})
deltapower_p = anova1(RdTot, allsubjects{2})
title('ANOVA: delta power')
xticklabels({'Normcog', 'MCI', 'PDD'})
gammapower_p = anova1(RgTot, allsubjects{2})
title('ANOVA: gamma power')
xticklabels({'Normcog', 'MCI', 'PDD'})
thetapower_p = anova1(RtTot, allsubjects{2})
title('ANOVA: theta power')
xticklabels({'Normcog', 'MCI', 'PDD'})

%% 5. Comparison of total power between bands and subjects

RaTot = sum(Ra, 2); %mean or sum? If you want to compare total power in the band, it should be sum
RbTot = sum(Rb, 2);
RdTot = sum(Rd, 2);
RgTot = sum(Rg, 2);
RtTot = sum(Rt, 2);

fiveBandTots = [RaTot, RbTot, RdTot, RgTot, RtTot];

[~, fiveBandTots_p, fiveBandTots_stats] = manova1(fiveBandTots, allsubjects{2}, 0.05)

figure, manovacluster(fiveBandTots_stats), title('MANOVA of within-band powers')

%% 6. Calculating Pillai's trace for MANOVA in case of very low p-values (<0.001)

%(Based on resource available at
%'https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Multivariate_Analysis_of_Variance-MANOVA.pdf')

Vs = sum((fiveBandTots_stats.eigenval)./(1+(fiveBandTots_stats.eigenval))) %Pillai's trace. Equivalent expression for it trace(fiveBandTots_stats.B*(fiveBandTots_stats.W + fiveBandTots_stats.B)^-1)

e = fiveBandTots_stats.dfW; %error term (within-class variability?)
h = fiveBandTots_stats.dfB; %hypothesis term (between-class variability?)

p = 5; %p is the total number of variables in the data
s = min(p, h);

m = (abs(p - h)-1)/2;
n = (e - p -1)/2;

F = (2*n + s + 1)*(Vs)/((2*m + s + 1)*(s - Vs)) %F-statistic approximation

1-fcdf(F, s*(2*m + s + 1), s*(2*m + s + 1)) %corresponding p-value...

%% 7. PCA and biplot

[pcaComps, pcaScores, ~, ~, pctExplained] = pca(fiveBandTots); %e contains the percent of variance explained by each eigenvector
pctExplained
biplot(pcaComps(:, 1:3), 'Scores',pcaScores(:, 1:3), 'VarLabels',{'comp1','comp2','comp3','comp4','comp5'}), axis vis3d

figure;
h = scatter(pcaScores(1:35, 1), pcaScores(1:35, 2), 16, 'cyan'); axis vis3d
hold on
scatter(pcaScores(36:72, 1), pcaScores(36:72, 2), 16, [1 .6 0])
scatter(pcaScores(72:end, 1), pcaScores(72:end, 2), 16, 'red')
hold off
legend('Normcog', 'MCI', 'PDD')
xlabel('PCA Component 1')
ylabel('PCA Component 2')
title('PCA analysis of MEG power-per-band')
