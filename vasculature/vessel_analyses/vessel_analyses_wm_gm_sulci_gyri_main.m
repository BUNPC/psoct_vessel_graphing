%% Vessel metrics - statistical analyses
%{
This script is for analyzing the length density in the dataset:
/projectnb/npbssmic/ns/Ann_Mckee_samples_10T
%}
clear; clc; close all;

%% Add top-level directory of code repository to path
% This allows Matlab to find the functions in the project folders

% Start in current directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Truncate path to reach top-level directory (psoct_vessel_graphing)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));

%% Data Directories and filenames
% Path to top-level directory
dpath = '/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/';

% Additional subfolder when applying minimum voxel count
% vox_min = 'vox_min_50';
vox_min = false;
% Volume and graph directories
voldir = '/dist_corrected/volume/ref/';
graphdir = ['/dist_corrected/volume/combined_segs/' ...
            'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];

% Update directory paths if using the minimum voxel
if vox_min
    graphdir = fullfile(graphdir, vox_min);
    mpath = fullfile(mpath, vox_min);
end

% Boolean for plotting violin plots in logarithmic scale
log_bool = true;

%% Import the metrics and heatmap structs
%%% Metrics - average value for each region
metrics = load(fullfile(mpath, 'metrics.mat'));
metrics = metrics.metrics;
subids = fields(metrics);
metrics_out = fullfile(mpath, 'metrics.mat');
% Average Metric Parameters. The tortuosity and diameter are continuous 
% distributions, which are tested with the LME model for the entire volume
params = {'length_density','branch_density','fraction_volume',...
            'tort_outliers','tortuosity','diameter'};
% Output filename for saving the p-value table
ptable_out = 'p_value_table.xls';

%%% Heatmaps struct - distribution of values for each region
roi_size = 1000;
hm_filename = append('heatmap_distro_',num2str(roi_size),'.mat');
hm = load(fullfile(mpath,hm_filename));
hm = hm.hm_distro;
% Average Heatmap Parameters
hm_params = {'vf','ld','bd','tort','diam'};
% Output filename for saving the heatmap p-value table
ptable_out_heatmap = append('p_value_table_heatmap_distro_',...
                            num2str(roi_size),'.xls');

%%% Brain regions for statistics
% Region of brain (excluding ratios)
regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
            'gm_gyri','wm_gyri'};
% Regions of brain (including ratios)
all_regions = {'tiss','gyri','sulci','gm','wm','gm_sulci','wm_sulci',...
               'gm_gyri','wm_gyri','sulci_gyri','wm_sulci_gyri',...
               'gm_sulci_gyri'};
% X-axis labels for violin plots and swarmcharts
xlabels = {'Volume','Gyri','Sulci','GM','WM','GM Sulci','WM Sulci',...
               'GM Gyri','WM Gyri'};

%%% Voxel dimensions (microns) and volume (cubic micron)
vox_dim = [12, 12, 15];
vox_vol = vox_dim(1) .* vox_dim(2) .* vox_dim(3);

%% Calculate N tortuosity outliers for each subject
% Calculate the outlier threshold from the healthy control group
% Then, calculate the number of outliers for each group (for each tissue
% region)

%%% Calculate outlier threshold based on the healthy control group
for ii = 1:length(regions)
    % Combine tortuosity measurements from all normal controls
    [~, ~, nc] = organize_metrics(metrics, subids,regions{ii}, 'tortuosity');
    
    % Calculate outlier cutoff = 3*MAD (scaled). This is taken from the
    % Mathworks website for calculating the threshold
    c = -1 / (sqrt(2)*erfcinv(3/2));
    mad3 = 3*c*mad(nc,1);
    cutoff = median(nc) + mad3;

    %%% Calculate number of outliers for each subject
    for j=1:length(subids)
        % Load tortuosity for each subject/region
        t = metrics.(subids{j}).(regions{ii}).tortuosity;
        n = sum(t > cutoff);
        % Save the number of outliers to struct
        metrics.(subids{j}).(regions{ii}).tort_outliers = n;
    end
end

%% Reorganize data and generate barcharts
% This script will reorganize the metrics struct by taking the region from
% each subject and combining them into a single substruct (e.g. taking the
% "tiss" region from each subject and creating a substruct named
% "metrics.tiss". This section can be skipped if it was already run.

%%% Barchart labels & titles
% Title of each bar chart
titles = {'Length Density','Branch Density','Volume Fraction',...
    'N Tortuosity Outliers'};
% Region titles
r_titles = {'Tissue', 'Gyri', 'Sulci', 'GM', 'WM',...
            'GM Sulci', 'WM Sulci', 'GM Gyri', 'WM Gyri'};
% Y-axis labels
ylabels = {'Length Density (mm/mm^2)','Branch Density (mm/mm^3)',...
            'Volume Fraction (a.u.)','Count'};
% Name of file to save
plot_names = {'length_density','branch_density','volume_fraction'...
    'tortuosity_outliers'};

%%% Iterate over each tissue region
for ii = 1:length(regions)
    %%% Organize the average metric for each subject
    % This includes the array of tortuosity values for each subject
    for j = 1:length(params)
        % Retrieve the average metric values for this region/parameter from
        % all AD, CTE, and NC subjects. Concatenate into arrays.
        [ad, cte, nc] = organize_metrics(metrics, subids,...
                                         regions{ii}, params{j});
        %%% Save parameter by the group to fascilitate statistical analyses
        metrics.(regions{ii}).(params{j}).ad = ad;
        metrics.(regions{ii}).(params{j}).cte = cte;
        metrics.(regions{ii}).(params{j}).nc = nc;
        % Create and save the bar chart (except for tortuosity
%         if j~=5
%             % Filename to save
%             plot_name = strcat(regions{ii},'_',plot_names{j});
%             % Title for bar chart
%             tstr = append(r_titles{ii},' - ', titles{j});
%             % Create bar chart and save output
%             bar_chart(ad, cte, nc, subids, mpath, tstr,'',...
%                       ylabels{j}, plot_name);
%             close;
%         end
    end

    %%% Organize the heatmap ROI values for each subject
    for j=1:length(hm_params)
        % Retrieve the heatmap ROI values for this region/parameter from
        % all AD, CTE, and NC subjects. Concatenate into arrays.
        [ad, cte, nc] = organize_metrics(hm,subids,...
                                         regions{ii},hm_params{j});
        % Save parameter by the group to fascilitate statistical analyses
        % in a later step.
        hm.(regions{ii}).(hm_params{j}).ad = ad;
        hm.(regions{ii}).(hm_params{j}).cte = cte;
        hm.(regions{ii}).(hm_params{j}).nc = nc;
    end
end

%% (Fig. 2) Violin Plots: Heatmap distribution intra-group differences

% Plot titles
vtitle_param = {'Volume Fraction','Length Density','Branch Density',...
    'Tortuosity','Diameter'};
vtitle_tissue = {'Entire Volume','Gyri','Sulci','GM','WM',...
                'GM Sulci','WM Sulci','GM Gyri','WM Gyri'};
% Y-axis labels
ylabels = {'Volume Fraction (a.u.)','Length Density (mm/mm^3)',...
            'Branch Density (1/mm^3)','Tortuosity (a.u.)',...
            'Diameter (\mum)'};

%%% Iterate over vascular metrics
for ii = 1:length(hm_params)
    %%% Iterate over brain regions
    for j = 1:length(regions)
        % Initialize array to store the values
        v_array = [];
        % Store disease/control group index (AD, CTE, HC)
        group_idx = {};
        % Store subject ID (AD_123,...)
        sub_idx = {};
        % Counters for number of violin plots for each group. Start with 1
        % to account for the whole group violin plot.
        n_ad = 1; n_cte = 1; n_hc = 1;
        
        %%% Iterate over subjects
        for k = 1:length(subids)
            % Retrieve array for subject
            sub_distro = hm.(subids{k}).(regions{j}).(hm_params{ii});
            % Create cell array for subject ID
            sub_str = repmat({subids{k}},[length(sub_distro),1]);
            % Create label cell array
            if contains(subids{k},'AD')
                sub_group = repmat({'AD'},[length(sub_distro),1]);
                n_ad = n_ad + 1;
            elseif contains(subids{k},'CTE')
                sub_group = repmat({'CTE'},[length(sub_distro),1]);
                n_cte = n_cte + 1;
            else
                sub_group = repmat({'HC'},[length(sub_distro),1]);
                n_hc = n_hc + 1;
            end
            % Concatenate the group and subject index arrays
            group_idx = vertcat(group_idx, sub_group);
            sub_idx = vertcat(sub_idx, sub_str);
            % Vertically concatenate the distro arrays into single vert array
            v_array = vertcat(v_array, sub_distro);
        end
        
        %%% Initialize table for violin plots
        vtable = table(group_idx,sub_idx, v_array);
    
        %%% Call function to generate violin plot
        % violin plot (no points)
        grpandplot(vtable,"v_array",yTitle=ylabels{ii},...
                   xFactor="sub_idx",...
                   cFactor="group_idx",...
                   showXLine=false,showVln=true,showBox=false,...
                   showPnt=false, showMean=false,log=log_bool,...
                   pntSize=5,w = 1,pntFillC='k',gap=0.6);
        title({vtitle_tissue{j},vtitle_param{ii}});
        set(gca, 'FontSize', 30)
        % Save output
        fname = append('heatmap_',hm_params{ii},'_',regions{j},'_violin');
        fout = fullfile(mpath,'heatmaps',fname);
        pause(0.1)
        saveas(gcf,fout,'png');
        close;
    end
end
%}

%% (Fig. 2) Violin Plots: Heatmap distribution inter-group differences

%%% Iterate over vascular metrics
for ii = 1:length(hm_params)
    % Initialize array to store the values
    v_array = [];
    % Store disease state index (AD, CTE, HC)
    group_idx = {};
    % Store brain region index ()
    region_idx = {};

    %%% Iterate over brain regions
    for j = 1:length(regions)
        % Vertically concatenate the [AD; CTE; NC] into single vert array
        ad = hm.(regions{j}).(hm_params{ii}).ad;
        cte = hm.(regions{j}).(hm_params{ii}).cte;
        nc = hm.(regions{j}).(hm_params{ii}).nc;
        v_array = vertcat(v_array, ad, cte, nc);
        n_samples = length(ad) + length(cte) + length(nc);

        % Initialize labels for comparisons within each group (AD, CTE, HC)
        ad = repmat({'AD'},[length(ad),1]);
        cte = repmat({'CTE'},[length(cte),1]);
        nc = repmat({'HC'},[length(nc),1]);
        group_idx = vertcat(group_idx, ad, cte, nc);

        % Initialize label for comparisons between brain regions
        region_idx = vertcat(region_idx,...
            repmat(cellstr(xlabels{j}),[n_samples,1]));
    end

    %%% Initialize table for violin plots
    vtable = table(region_idx, group_idx, v_array);

    %%% Call function to generate violin plot
    % Group each AD/HC/CTE comparison by tissue
    grpandplot(vtable,"v_array", yTitle = ylabels{ii},...
        xFactor="region_idx", cFactor = "group_idx", xOrder = xlabels,...
        showXLine = true, showVln = true, showBox = false, showMean=false,...
        showPnt = false, showNum = false, numYPos = 500, pntSize=5,...
        gap=0.6, log=log_bool);
    set(gca, 'FontSize', 30)
    % Save output
    fout = fullfile(mpath,'/heatmaps/',append('heatmap_',hm_params{ii},'_violin'));
    pause(0.1)
    saveas(gcf,fout,'png');
    close;
end
%}

%% Fig.2 - Heatmap Statistical Hypothesis Testing (LME model)
% Significant Difference threshold
alpha = 0.05;
% Calculate stats
hm_stats = calc_heatmap_stats(hm, regions, hm_params, subids, alpha,...
                             mpath, ptable_out_heatmap);

%% Calculate ratio of sulci/gyri
% This cannot be calculated for the tortuosity since it's an array of
% values.

% Cell array of parameters to compute ratio 
ratio_params = {'length_density','branch_density','fraction_volume',...
    'tort_outliers'};
ratio_params = params;
groups = {'ad','cte','nc'};

% Iterate over each parameter
for ii = 1:length(ratio_params)
    % Iterate over groups
    for j = 1:length(groups)
        %%% Retrieve the parameter for each region
        sulci = metrics.sulci.(ratio_params{ii}).(groups{j});
        gyri = metrics.gyri.(ratio_params{ii}).(groups{j});
        wm_sulci = metrics.wm_sulci.(ratio_params{ii}).(groups{j});
        gm_sulci = metrics.gm_sulci.(ratio_params{ii}).(groups{j});
        wm_gyri = metrics.wm_gyri.(ratio_params{ii}).(groups{j});
        gm_gyri = metrics.gm_gyri.(ratio_params{ii}).(groups{j});
        
        %%% Compute the ratios
        % sulci / gyri ratios
        sulci_gyri = sulci ./ gyri;
        % WM sulci / gyri ratios
        wm_sulci_gyri = wm_sulci ./ wm_gyri;
        % GM sulci / gyri ratios
        gm_sulci_gyri = gm_sulci ./ gm_gyri;

        %%% Save to the struct
        metrics.sulci_gyri.(ratio_params{ii}).(groups{j}) = sulci_gyri;
        metrics.wm_sulci_gyri.(ratio_params{ii}).(groups{j}) = wm_sulci_gyri;
        metrics.gm_sulci_gyri.(ratio_params{ii}).(groups{j}) = gm_sulci_gyri;
    end
end

%% Fig. 3: Statistical Hypothesis Testing (average for each region)

% Trend threshold
trend = 0.10;
% Significant Difference threshold
alpha = 0.05;
% Calculate stats
pstats = calc_avg_stats(metrics, all_regions, params, alpha,...
                        trend, mpath, ptable_out);
% Save statistics
stats_fout = fullfile(mpath, 'stats.mat');
save(stats_fout, 'pstats','-v7.3');

%% Fig. 3 - Swarmchart for average values (inter-group)

% y-axis labels
ylabels = {'Length Density (mm/mm^3)','Branch Density (1/mm^3)',...
            'Volume Fraction (a.u.)','N Outliers',...
            'Tortuosity (a.u.)','Diameter (\mum)'};

for ii = 1:length(params)
    % Initialize array to store the values
    v_array = [];
    % Store disease state index (AD, CTE, HC)
    group_idx = {};
    % Store brain region index ()
    region_idx = {};

    %%% Iterate over brain regions
    for j = 1:length(regions)
        % Vertically concatenate the [AD; CTE; NC] into single vert array
        ad = metrics.(regions{j}).(params{ii}).ad;
        cte = metrics.(regions{j}).(params{ii}).cte;
        nc = metrics.(regions{j}).(params{ii}).nc;
        v_array = vertcat(v_array, ad, cte, nc);
        n_samples = length(ad) + length(cte) + length(nc);

        % Initialize labels for comparisons within each group (AD, CTE, HC)
        ad = repmat({'AD'},[length(ad),1]);
        cte = repmat({'CTE'},[length(cte),1]);
        nc = repmat({'HC'},[length(nc),1]);
        group_idx = vertcat(group_idx, ad, cte, nc);

        % Initialize label for comparisons between brain regions
        region_idx = vertcat(region_idx,...
            repmat(cellstr(xlabels{j}),[n_samples,1]));
    end

    %%% Initialize table for violin plots
    vtable = table(region_idx, group_idx, v_array);

    %%% Call function to generate violin plot
    % Group each AD/HC/CTE comparison by tissue
    grpandplot(vtable,"v_array", yTitle=ylabels{ii},...
        xFactor="region_idx", cFactor="group_idx", xOrder=xlabels,...
        showXLine=true, showVln=false, showBox=false, showMean=true,...
        showPnt=true, pntOnTop=true,jitter=false,...
        showNum=false, numYPos=500, pntSize=40,...
        gap=0.6, log=log_bool);
    set(gca, 'FontSize', 30)
    % Save output
    if ~exist(fullfile(mpath,'/swarmchart/'))
        mkdir(fullfile(mpath,'/swarmchart/'));
    end
    fout = fullfile(mpath,'/swarmchart/',...
                    append('swarmchart_',params{ii}));
    pause(0.1);
    saveas(gcf,fout,'png');
    close;
end

%% Fig. 3 Generate boxplots for average values (deprecated)
%{
% Hui decided she no longer wanted boxplots, so this code is deprecated.
% Create boxplots for the metrics that have significance:
% Length Density: WM sulci (AD<HC, CTE<HC), sulci (AD<HC)
% Branch Density: WM sulci (AD<HC)

% Groups for binning the violin plots
groups = {'AD','CTE','HC'};
% y-axis labels
ylabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)'};
% Box/Whisker Plot Titles
bp_region_title = {'White Matter Sulci','Sulci'};
bp_param_title = {'Length Density','Branch Density'};

% Iterate over parameters
for ii = 1:length(bp_params)
    % Initialize array to store the values
    v_array = [];
    % Store disease state index (AD, CTE, HC)
    group_idx = {};
    % Store brain region index ()
    region_idx = {};
    % Iterate over regions
    for j = 1:length(bp_regions)
        % Vertically concatenate the [AD; CTE; NC] into single vert array
        ad = metrics.(bp_regions{j}).(bp_params{ii}).ad;
        cte = metrics.(bp_regions{j}).(bp_params{ii}).cte;
        nc = metrics.(bp_regions{j}).(bp_params{ii}).nc;
        v_array = vertcat(ad, cte, nc);
        n_samples = length(ad) + length(cte) + length(nc);
        
        % Initialize labels for comparisons within each group (AD, CTE, HC)
        ad = repmat({'AD'},[length(ad),1]);
        cte = repmat({'CTE'},[length(cte),1]);
        nc = repmat({'HC'},[length(nc),1]);
        group_idx = vertcat(ad, cte, nc);
        
        %%% Initialize table for violin plots
        vtable = table(group_idx, v_array);
    
        %%% Call function to generate violin plot
        % Group each AD/HC/CTE comparison by tissue
        grpandplot(vtable,"v_array", yTitle = ylabels{ii},...
                    xFactor = 'group_idx',...
                    cFactor = 'group_idx',...
                    pntFillC = 'k',...
                    showBox = true, showPnt = true, pntSize=40, w=1);
        set(gca, 'FontSize', 50)
        % Create title
        bp_title = {bp_region_title{j},bp_param_title{ii}};
        title(bp_title);
        % Save output
        fout = fullfile(mpath,'/box_whisker/',...
                append(bp_params{ii},'_',bp_regions{j},'_box_whisker'));
        saveas(gcf,fout,'png');
        close;
    end
end
%}

%% Violin plots of average values for each sample (deprecated)
% This section will iterate over each vascular metric and generate a 1xm
% cell array, where each column contains the vascular metric values for a
% particular brain region. Then, this section will generate a violin plot
% for each vascular metric, where the x-axis will represent different brain
% regions. For each brain region, the AD, CTE, NC will be grouped together

%{
% Y-axis labels
ylabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)',...
            'Volume Fraction (a.u.)','Tortuosity (a.u.)'};
% X-axis labels for groups
xlabels = {'Volume','Gyri','Sulci','GM','WM','GM Sulci','WM Sulci',...
               'GM Gyri','WM Gyri'};
% Regions of the brain for calculating ratios
ratio_regions = {'sulci_gyri','wm_sulci_gyri','gm_sulci_gyri'};
ratio_ylabels = {'Length Density (\mum^-^2)','Branch Density (\mum^-^3)',...
            'Volume Fraction (a.u.)'};
ratio_xlabels = {'Sulci / Gyri', 'WM Sulci / Gyri', 'GM Sulci / Gyri'};

%%% Iterate over vascular metrics
for ii = 1:length(params)   
    % Initialize array to store the values
    v_array = [];
    % Store disease state index (AD, CTE, HC)
    group_idx = {};
    % Store brain region index ()
    region_idx = {};

    %% Iterate over brain regions (exclude ratios)
    for j = 1:length(regions)
        % Vertically concatenate the [AD; CTE; NC] into single vert array
        ad = metrics.(regions{j}).(params{ii}).ad;
        cte = metrics.(regions{j}).(params{ii}).cte;
        nc = metrics.(regions{j}).(params{ii}).nc;
        v_array = vertcat(v_array, ad, cte, nc);
        n_samples = length(ad) + length(cte) + length(nc);

        % Initialize labels for comparisons within each group (AD, CTE, HC)
        ad = repmat({'AD'},[length(ad),1]);
        cte = repmat({'CTE'},[length(cte),1]);
        nc = repmat({'HC'},[length(nc),1]);
        group_idx = vertcat(group_idx, ad, cte, nc);

        % Initialize label for comparisons between brain regions
        region_idx = vertcat(region_idx,...
            repmat(cellstr(xlabels{j}),[n_samples,1]));
    end

    %%% Initialize table for violin plots
    vtable = table(region_idx, group_idx, v_array);

    %%% Call function to generate violin plot
    % Group each AD/HC/CTE comparison by tissue
    grpandplot(vtable,"v_array", yTitle = ylabels{ii},...
        xFactor="region_idx", cFactor = "group_idx", xOrder = xlabels,...
        showXLine = true, showVln = true, showBox = false, showMean=true,...
        showPnt = false, showNum = false, numYPos = 500, pntSize=5);
    set(gca, 'FontSize', 30)
    % Save output
    fout = fullfile(mpath,append(params{ii},'_violin'));
    saveas(gcf,fout,'png');
    close;
    %% Iterate over brain regions (only ratios)
    if strcmp(params{ii}, 'tortuosity')
        continue
    else
        % Initialize array to store the values
        v_array = [];
        % Store disease state index (AD, CTE, HC)
        group_idx = {};
        % Store brain region index ()
        region_idx = {};
        
        % Iterate over ratio regions
        for j = 1:length(ratio_regions)
            % Vertically concatenate [AD; CTE; NC] into single vert array
            ad = metrics.(ratio_regions{j}).(params{ii}).ad;
            cte = metrics.(ratio_regions{j}).(params{ii}).cte;
            nc = metrics.(ratio_regions{j}).(params{ii}).nc;
            v_array = vertcat(v_array, ad, cte, nc);
            n_samples = length(ad) + length(cte) + length(nc);
    
            % Initialize labels for comparisons within group (AD, CTE, HC)
            ad = repmat({'AD'},[length(ad),1]);
            cte = repmat({'CTE'},[length(cte),1]);
            nc = repmat({'HC'},[length(nc),1]);
            group_idx = vertcat(group_idx, ad, cte, nc);
    
            % Initialize label for comparisons between brain regions
            region_idx = vertcat(region_idx,...
                repmat(cellstr(ratio_xlabels{j}),[n_samples,1]));
        end
    
        %%% Initialize table for violin plots
        vtable = table(region_idx, group_idx, v_array);
    
        %%% Call function to generate violin plot
        % Group each AD/HC/CTE comparison by tissue
        grpandplot(vtable,"v_array", yTitle=ratio_ylabels{ii},...
            xFactor="region_idx", cFactor="group_idx",...
            xOrder = ratio_xlabels, showXLine = true, showVln = true,...
            showBox = false, showMean=true,...
            showPnt = false, showNum = false, numYPos = 500, pntSize=5);
        set(gca, 'FontSize', 30)
        % Save output
        fout = fullfile(mpath,append(params{ii},'_ratio_violin'));
        saveas(gcf,fout,'png');
        close;
    end

end
%}