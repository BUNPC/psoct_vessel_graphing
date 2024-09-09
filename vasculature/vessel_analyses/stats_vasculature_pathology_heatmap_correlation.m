%% Measure correlation b/w pathology and vascular heatmap matrices
% Overview:
%   - overlay the vascular heatmap with the registered pathology
%   - iterate over the vascular heatmap ROIs (WM and GM separately)
%       - measure the average pathology value for each vascular ROI
%       - add [x,y] pair to matrix
%   - Statistics:
%       - spearman's rho corr. [rho, p] = corr(X,Y,'Type','Spearman')

%% Add top-level directory of code repository to path
clear; clc; close all;
% Print current working directory
mydir  = pwd;
% Find indices of slashes separating directories
if ispc
    idcs = strfind(mydir,'\');
elseif isunix
    idcs = strfind(mydir,'/');
end
% Remove the two sub folders to reach parent
% (psoct_human_brain\vasculature\vesSegment)
topdir = mydir(1:idcs(end-1));
addpath(genpath(topdir));
% Set maximum number of threads equal to number of threads for script
NSLOTS = str2num(getenv('NSLOTS'));
maxNumCompThreads(NSLOTS);

%% Initialize directories, filenames, parameters
%%% Directories 
% Metrics output path
mpath = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/'];
% Registered staining heatmaps directory
path_reg = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
    'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/heatmaps/' ...
    'heatmaps_pathology_registration/'];

%%% Vasculature heatmap matrix filename
% ROI cube side (microns)
cube_side = 1000;
% Load according to size of ROI cube
hm_fname = append('heatmap_ab_ptau_',num2str(cube_side),'.mat');
% Load the vascular heatmap matrix
hm = load(fullfile(mpath,hm_fname));
hm = hm.heatmap;
subid = fields(hm);

%%% Initialize filenames of registered pathology / heatmaps.
% Some of the subjects had the patholoy registered with the automated
% "demons" algorithm and other were registered with the Maritos manual
% landmark software. This section will specify this delineation.
%
% Subjects with heatmaps in TIF
ftif_ab_path = 'Ab_registered.tif';
ftif_ab_mask = 'Ab_mask_registered.tif';
ftif_pt_path = 'AT8_registered.tif';
ftif_pt_mask = 'AT8_mask_registered.tif';
subs_tif = {'AD_20832','AD_20969','CTE_6489','CTE_6912',...
    'NC_6839','NC_21499'};
% Subjects with heatmaps in .MAT (the actual filenames have the subject ID
% appended as a prefix (i.e. [sub_id]_[file_name])
fmat_ab_path = 'Ab_path_registered.mat';
fmat_ab_mask = 'Ab_path_mask_registered.mat';
fmat_pt_path = 'AT8_path_registered.mat';
fmat_pt_mask = 'AT8_path_mask_registered.mat';
subs_mat = {'AD_10382','AD_21354','AD_21424','CTE_7019','CTE_7126',...
    'NC_8095'};

%%% Subvolume parameters
% Isotropic cube length (microns)
cube_side = 1000;
% Size of each OCT voxel (microns)
vox = [12, 12, 15];
% Whether to plot non-normalized heatmaps for each depth
viz_individual = false;

%%% Compute number of voxels in x,y,z for each cube
n_x = floor(cube_side ./ vox(1));
n_y = floor(cube_side ./ vox(2));
n_z = floor(cube_side ./ vox(3));
% Calculate the size of each cube (in voxels)
cube_vol_vox = n_x * n_y * n_z;
% Calculate the size of each cube (in cubic microns)
cube_vol_um = cube_vol_vox * vox(1) * vox(2) * vox(3);

%%% Index for structs
% Pathology
pidx = {'ab','pt'};
% masks (gray matter and white matter)
midx = {'gm','wm'};

%%% Struct for storing pairs of values (heatmap, pathology)
pairs = struct();

%% Spearman's Rho Correlation Coefficient b/w pathology & vasculature
% Load the vascular heatmap and pathology heatmaps, recreate the ROIs from
% the vascular heatmaps, generate an x-y pair for each ROI (x = vascular
% ROI, y = pathology ROI), perform Spearman's on the matrix pairs for each
% subject and stain.

% Iterate over each subject
for ii = 1:length(subid)
    %% Load the registered pathology and retrieve vascular heatmaps
    % Retrieve the current subject ID
    sub = subid{ii};
    % Determine which registration method was used on this subject 
    if any(ismember(subs_tif,sub))
        % TIF: Load the A-beta and p-tau pathology heatmaps
        ab = TIFF2MAT(fullfile(path_reg,sub,ftif_ab_path));
        ab_mask = TIFF2MAT(fullfile(path_reg,sub,ftif_ab_mask));
        pt = TIFF2MAT(fullfile(path_reg,sub,ftif_pt_path));
        pt_mask = TIFF2MAT(fullfile(path_reg,sub,ftif_pt_mask));
        % Convert from RGB to grayscale (all depths are equivalent
        ab = ab(:,:,1); ab_mask = logical(ab_mask(:,:,1));
        pt = pt(:,:,1); pt_mask = logical(pt_mask(:,:,1));
    else
        % MAT: Load the A-beta and p-tau pathology heatmaps
        ab = load(fullfile(path_reg,sub,append(sub,'_',fmat_ab_path)));
        ab_mask = load(fullfile(path_reg,sub,append(sub,'_',fmat_ab_mask)));
        pt = load(fullfile(path_reg,sub,append(sub,'_',fmat_pt_path)));
        pt_mask = load(fullfile(path_reg,sub,append(sub,'_',fmat_pt_mask)));
        % Load the fields of each struct
        ab = ab.path_registered; pt = pt.path_registered;
        ab_mask = logical(ab_mask.path_mask_registered);
        pt_mask = logical(pt_mask.path_mask_registered);
    end

    % Retrieve the vascular heatmaps/masks from the "hm" struct
    vasc = hm.(sub);
    masks = struct();
    masks.gm = logical(vasc.mask_gm);
    masks.wm = logical(vasc.mask_wm);
    hm_vf = vasc.vf;
    hm_bd = vasc.bd;
    hm_ld = vasc.ld;
    hm_tr = vasc.tort;
    hm_dm = vasc.diam;

    % Combine the masks and pathology into struct
    path = struct();
    path.ab.mask = ab_mask;
    path.pt.mask = pt_mask;
    path.ab.hm = ab;
    path.pt.hm = pt;

    
    %% Iterate over the pathology stains
    % Iteration 1 = A-beta, Iteration 2 = p-tau
    for j=1:2
        %%% Iterate over the white matter and gray matter
        % 1 = GM, 2 = WM
        for k=1:2
            % Initialize matrices for [x,y] pairs for each metric
            n_pair_x = size(1:n_x:size(hm_vf,1),2);
            n_pair_y = size(1:n_y:size(hm_vf,2),2);
            n_pairs = n_pair_x .* n_pair_y;
            vf_pairs = zeros(n_pairs,2);
            ld_pairs = zeros(n_pairs,2);
            bd_pairs = zeros(n_pairs,2);
            tr_pairs = zeros(n_pairs,2);
            dm_pairs = zeros(n_pairs,2);
            % Track the pair index in the nested for loops
            idx = 1;
            % Combine the pathology / vascular masks
            mask = masks.(midx{k});
            mask = mask(:,:,j) .* path.(pidx{j}).mask;
            % Retrieve the respective pathology heatmap matrix
            hm_path = path.(pidx{j}).hm;
            
            %% Iterate over ROIs
            % Iterate over rows
            for x = 1:n_x:size(hm_vf,1)
                % Iterate over columns
                for y = 1:n_y:size(hm_vf,2)
                    %% Set the x/y bounds 
                    % Initialize end indices for each axis
                    xf = x + n_x - 1;
                    yf = y + n_y - 1;
                    % Take minimum of matrix dimensions and end indices
                    xf = min(xf, size(hm_vf,1));
                    yf = min(yf, size(hm_vf,2));
                    
                    %% Determine if combined mask contain tissue
                    if any(mask((x:xf),(y:yf)))
                        % Take the minimum vascular heatmap value. This should
                        % be uniform across the ROI, so minimum is arbitrary
                        vf_tmp = min(hm_vf((x:xf), (y:yf),j),[],"all");
                        bd_tmp = min(hm_bd((x:xf), (y:yf),j),[],"all");
                        ld_tmp = min(hm_ld((x:xf), (y:yf),j),[],"all");
                        tr_tmp = min(hm_tr((x:xf), (y:yf),j),[],"all");
                        dm_tmp = min(hm_dm((x:xf), (y:yf),j),[],"all");
                        
                        % Take average pathology value within ROI
                        path_tmp = mean(hm_path((x:xf), (y:yf)),'all');
    
                        % Create x-y pair for each metric
                        vf_pairs(idx,:) = [vf_tmp, path_tmp];
                        bd_pairs(idx,:) = [bd_tmp, path_tmp];
                        ld_pairs(idx,:) = [ld_tmp, path_tmp];
                        tr_pairs(idx,:) = [tr_tmp, path_tmp];
                        dm_pairs(idx,:) = [dm_tmp, path_tmp];
                        
                        % Iterate counter index
                        idx = idx + 1;
                    else
                        continue
                    end
                end
            end
            % Retain pairs up to the index
            pairs.(sub).(midx{k}).(pidx{j}).vf = vf_pairs(1:idx,:);
            pairs.(sub).(midx{k}).(pidx{j}).bd = bd_pairs(1:idx,:);
            pairs.(sub).(midx{k}).(pidx{j}).ld = ld_pairs(1:idx,:);
            pairs.(sub).(midx{k}).(pidx{j}).tr = tr_pairs(1:idx,:);
            pairs.(sub).(midx{k}).(pidx{j}).dm = dm_pairs(1:idx,:);
        end
    end
    sprintf('Finished subject %s\n',sub)
end

%% Statistics (spearman's rho correlations)
% Initialize row labels for table of spearman's rho
Stat = {'rho';'p-value'};

% Initialize struct to store the spearman's rho
spear = struct();

% Iterate over each subject
for ii = 1:length(subid)
    % Retrieve subject 
    sub = subid{ii};
    % Iterate over the stain
    for j=1:2
        % Iterate over GM, WM
        for k=1:2
            % Load the pairs
            vf = pairs.(sub).(midx{k}).(pidx{j}).vf;
            bd = pairs.(sub).(midx{k}).(pidx{j}).bd;
            ld = pairs.(sub).(midx{k}).(pidx{j}).ld;
            tr = pairs.(sub).(midx{k}).(pidx{j}).tr;
            dm = pairs.(sub).(midx{k}).(pidx{j}).dm;
            % Call Spearman's
            [vf_rho, vf_p] = corr(vf(:,1),vf(:,2),'type','Spearman');
            [bd_rho, bd_p] = corr(bd(:,1),bd(:,2),'type','Spearman');
            [ld_rho, ld_p] = corr(ld(:,1),ld(:,2),'type','Spearman');
            [tr_rho, tr_p] = corr(tr(:,1),tr(:,2),'type','Spearman');
            [dm_rho, dm_p] = corr(dm(:,1),dm(:,2),'type','Spearman');
    
            %%% Assign to struct
            % Volume fraction
            spear.(sub).(midx{k}).(pidx{j}).vf.rho = vf_rho;
            spear.(sub).(midx{k}).(pidx{j}).vf.p = vf_p;
            % Length density
            spear.(sub).(midx{k}).(pidx{j}).ld.rho = ld_rho;
            spear.(sub).(midx{k}).(pidx{j}).ld.p = ld_p;
            % Branch density
            spear.(sub).(midx{k}).(pidx{j}).bd.rho = bd_rho;
            spear.(sub).(midx{k}).(pidx{j}).bd.p = bd_p;
            % Tortuosity
            spear.(sub).(midx{k}).(pidx{j}).tr.rho = tr_rho;
            spear.(sub).(midx{k}).(pidx{j}).tr.p = tr_p;
            % Diameter
            spear.(sub).(midx{k}).(pidx{j}).dm.rho = dm_rho;
            spear.(sub).(midx{k}).(pidx{j}).dm.p = dm_p;
    
    
            % Print to console if any of the p values are below 0.05
            if vf_p < 0.05
                sprintf('SIG: subject %s, path %s, met VF, p=%f, rho=%f',...
                    sub,pidx{j},vf_p,vf_rho)
            elseif ld_p < 0.05
                sprintf('SIG: subject %s, path %s, met LD, p=%f, rho=%f',...
                    sub,pidx{j},ld_p,ld_rho)
            elseif bd_p < 0.05
                sprintf('SIG: subject %s, path %s, met BD, p=%f, rho=%f',...
                    sub,pidx{j},bd_p,ld_rho)
            end
            %%% Add spearman's to table
            VolumeFraction = [vf_rho; vf_p];
            LengthDensity = [ld_rho; ld_p];
            BranchDensity = [bd_rho; bd_p];
            Tortuosity = [tr_rho; tr_p];
            Diameter = [dm_rho; dm_p];
            ptable = make_heatmap_ptable(Stat,VolumeFraction,LengthDensity,...
                                        BranchDensity,Tortuosity,Diameter);
            %%% Save to table
            % Name of sheet
            sheet = append(sub,'_',midx{k},'_',pidx{j});
            % Create output filename
            table_out = fullfile(mpath, 'p_value_spearmans.xls');
            writetable(ptable, table_out, 'Sheet',sheet);
        end
    end
end

%% Create single table of all subjects / pathologies / metrics
% Count number of subjects, metrics, and pathologies
nsub = length(fields(spear));
nmet = length(fields(spear.AD_10382.gm.ab));
npath = length(pidx);
% Initialize matrix to store the p-values and rho values
gm_rho_p = zeros((nsub .* npath),nmet);
wm_rho_p = zeros((nsub .* npath),nmet);

% Iterate over GM, WM
for k=1:2
    % Index to track row
    idx = 1;
    % Iterate over each subject
    for ii = 1:length(subid)
        % Retrieve subject 
        sub = subid{ii};
        % Iterate over the stain
        for j=1:2
            % Retrieve the rho and p-value for each metric
            vf = spear.(subid{ii}).(midx{k}).(pidx{j}).vf;
            ld = spear.(subid{ii}).(midx{k}).(pidx{j}).ld;
            bd = spear.(subid{ii}).(midx{k}).(pidx{j}).bd;
            tr = spear.(subid{ii}).(midx{k}).(pidx{j}).tr;
            dm = spear.(subid{ii}).(midx{k}).(pidx{j}).dm;
            % Add to matrix
            if strcmp(midx{k},'gm')
                gm_rho_p(idx,:) = [vf.rho, ld.rho, bd.rho, tr.rho, dm.rho];
                gm_rho_p(idx+1,:) = [vf.p, ld.p, bd.p, tr.p, dm.p];
            else
                wm_rho_p(idx,:) = [vf.rho, ld.rho, bd.rho, tr.rho, dm.rho];
                wm_rho_p(idx+1,:) = [vf.p, ld.p, bd.p, tr.p, dm.p];
            end
            % Iterate row index
            idx = idx + 2;
        end
    end
end

%%% Create table
% Strings for the pathology
path_cell = repmat({'A-beta';'A-beta';'p-tau';'p-tau'},[nsub,1]);
% Strings for rho or p-value
rho_p_cell = repmat({'rho';'p-value'},[nsub*npath,1]);
% Strings for the group (AD, CTE, HC)
ad_cell = repmat({'AD'},[5.*4,1]);
cte_cell = repmat({'CTE'},[4.*4,1]);
hc_cell = repmat({'HC'},[3.*4,1]);
sub_cell = vertcat(ad_cell, cte_cell, hc_cell);
% Combine GM rho and p-values into table
vf = gm_rho_p(:,1);
ld = gm_rho_p(:,2);
bd = gm_rho_p(:,3);
tr = gm_rho_p(:,4);
dm = gm_rho_p(:,5);
gm_table = table(sub_cell, path_cell, rho_p_cell,vf,ld,bd,tr,dm);
% Combine WM rho and p-values into table
vf = wm_rho_p(:,1);
ld = wm_rho_p(:,2);
bd = wm_rho_p(:,3);
tr = wm_rho_p(:,4);
dm = wm_rho_p(:,5);
wm_table = table(sub_cell, path_cell, rho_p_cell,vf,ld,bd,tr,dm);
% Write to spreadsheet
table_out = fullfile(mpath, 'p_value_spearmans.xls');
writetable(gm_table, table_out, 'Sheet', 'gm_combined');
writetable(wm_table, table_out, 'Sheet', 'wm_combined');

%% Fig. 4: Scatter plot of pathology vs. vasculature
% The struct "pairs" contains a matrix of [x,y] pairs for each subject,
% pathology, and vascular metric (i.e. pairs.[sub].[ab/pt].[vf/bd/ld])
% The first column is the vascular metrix (x-axis) within the ROI, and the
% second column is the pathology (y-axis) within the same ROI.

% Vascular Metrics
vmets = {'vf','bd','ld','tr','dm'};

% X-label cell array
xlabels = {'Volume Fraction (a.u.)','Branch Density','Length Density',...
    'Tortuosity','Diameter'};
ylabels = {'[A-beta] (a.u.)','[p-tau] (a.u.)'};

% Iterate over subject IDs
for ii = 1:length(subid)
    % Iterate over pathologies
    for j=1:length(pidx)
        % Iterate over vascular metrics
        for k=1:length(vmets)
            % Iterate over gray matter and white matter masks
            for m=1:2
                % Extract the pair values
                pair = pairs.(subid{ii}).(midx{m}).(pidx{j}).(vmets{k});
                % Set the x-axis label (vascular metrics)
                xl = xlabels{k};
                % Set the y-axis label (pathology)
                yl = ylabels{j};
                % Set title string
                tstr = append(subid{ii},' ',midx{m},' ',pidx{j},' ',vmets{k});
                % Set the filename
                fname = append(subid{ii},'_',midx{m},'_',pidx{j},'_',vmets{k},'_scatter.png');
                % Call plotting function
                lin_reg_plot(pair,subid{ii},tstr,xl,yl,path_reg,fname)
            end
        end
    end
end

%% Fig. 4: Generate figure of registered pathology heatmap 
% Figure for the manuscript of the AD/CTE/HC vasculature study
% A subset of samples were registered with the automated demons method, and
% the output was a .MAT file. The rest required a manual landmark method to
% register, and the output was a .TIF. This section will handle each case.
%
% Automated (.MAT): AD_10382, AD_21354, AD_21424, CTE_7019, CTE_7126
% Manual (.TIF): AD_20832, AD_20969, CTE_6489, CTE_6912, NC_8095, NC_6839,
%               NC_21499
% TODO:
%   - rotate the heatmaps that are at an angle so that gyri are on top

%%% Initialize filenames (excluding AB or AT8 prefix)
% Data path to registered mask and pathology
reg_path = ['/projectnb/npbssmic/ns/Ann_Mckee_samples_55T/metrics/' ...
        'gsigma_1-3-5_2-3-4_3-5-7_5-7-9_7-9-11/p18/heatmaps/' ...
        'heatmaps_pathology_registration'];
% Filename prefix for amyloid beta and AT8
pidx = {'Ab','AT8'};
% Subject list registered with automated demons
automated = {'AD_10382', 'AD_21354', 'AD_21424', 'CTE_7019', 'CTE_7126',...
             'NC_8095'};

%%% Load heatmaps and masks and make figures
% Iterate over subject IDs
for ii = 1:length(subid)
    % Iterate over pathologies
    for j=1:length(pidx)
        %%% Load file, depending on registration method
        % Check if automated registration method used (.MAT)
        if any(contains(automated, subid{ii}))
            % Define filenames
            mask = append(subid{ii},'_',pidx{j},'_path_mask_registered.mat');
            path = append(subid{ii},'_',pidx{j},'_path_registered.mat');
            % Load the registered pathology tissue mask
            mask = load(fullfile(reg_path,subid{ii},mask));
            mask = mask.path_mask_registered;
            % Load the registered pathology heatmap
            path = load(fullfile(reg_path,subid{ii},path));
            path = path.path_registered;
        % Manual landmark method (.TIF)
        else
            % Define filenames
            mask = append(pidx{j},'_mask_registered.tif');
            path = append(pidx{j},'_registered.tif');
            % Load the registered pathology tissue mask
            mask = TIFF2MAT(fullfile(reg_path,subid{ii},mask));
            path = TIFF2MAT(fullfile(reg_path,subid{ii},path));
            % Take one depth from stack
            mask = mask(:,:,1);
            path = path(:,:,1);
        end

        %%% Initialize the title and output filenames
        % Output path
        fpath = fullfile(reg_path,subid{ii});
        % Figure Name
        figname = append(subid{ii},'_',pidx{j},'_path_registered.png');
        title_str = append(subid{ii},' ',pidx{j});
        % Plot heatmap
        plot_save_heatmap(1,path,0,[0,256],mask,title_str,'(a.u.)',...
                         fpath,figname)
    end
end

%% Plot and save the heat maps
function plot_save_heatmap(Ndepths, heatmaps, flip_cbar, colorbar_range,...
    masks, tstr, cbar_label, dpath, fname)
% PLOT_SAVE_HEATMAP use imagesc and set background = 0
% INPUT
%   Ndepths (int): number of depths in z dimension
%   heatmaps (double matrix): heatmaps of vascular metric
%   flip_cbar (logical): reverse the direction of the colorbar
%   colorbar_range (double array): [min, max]
%   masks (double): tissue mask (1=tissue, 0=other)
%   tstr (string): figure title
%   cbar_label (string): colorbar label
%   dpath (string): data dicrectory path
%   fname (string): name of figure to save

%%% Set the number of depths to iterate for each heatmap
% If the number of depths is not specified, then set it equal to the number
% of z dimensions.
if isempty(Ndepths)
    Ndepths = size(heatmaps,3);
end
% Set fontsize for the heatmap figure
fontsize = 40;

%%% Iterate over frames in z dimension
for d = 1:Ndepths
    %%% Heatmap of j_th frame from the length density
    fh = figure();
    fh.WindowState = 'maximized';
    % If there are multiple heatmaps in the matrix
    if size(heatmaps,3) > 1
        heatmap = heatmaps(:,:,d);
    % Here it is just a single frame of a heatmap
    else
        heatmap = heatmaps;
    end
    % Initialize heatmap
    h = imagesc(heatmap);

    %%% Initialize colorbar
    % If the colorbar_range is passed in, then extract min & max
    if ~isempty(colorbar_range)
        cmap_min = colorbar_range(1);
        cmap_max = colorbar_range(2);
    % Otherwise, set limits from the current heatmap
    else
        % Min = Lowest value also greater than zero
        cmap_min = min(heatmap(heatmap(:)>0));
        % Max = Find the 95th percentile for upper limit
        cmap_max = prctile(heatmap(:), 95);
    end
    % Initialize the colormap limits
    cmap = jet(256);
    clim(gca, [cmap_min, cmap_max]);
    % Initialize colormap and colorbar
    if flip_cbar
        colormap(flipud(cmap));
    else
        colormap(cmap);
    end
    c = colorbar;

    %%% Apply tissue mask to the heatmap to remove background
    alpha_mask = double(masks(:,:,d));
    set(h, 'AlphaData', alpha_mask);

    %%% Configure figure parameters
    % Update title string with specific pathology
    pathology = {'A-Beta','p-tau'};
    if size(heatmaps,3) > 1
        title_str = append(tstr, ' ', pathology{d});
    else
        title_str = tstr;
    end
    title(title_str,'Interpreter','none');
    set(gca, 'FontSize', fontsize);
    % Label the colorbar    
    c.Label.String = cbar_label;
    % Offset colorbar label to the right of colorbar
    c.Label.Position = [10 (cmap_max - (cmap_max-cmap_min)/2)];
    c.Label.Rotation = 270;
    % Increase fontsize of colorbar
    c.FontSize = 40;
    % Remove x and y tick labels
    set(gca,'Yticklabel',[]);
    set(gca,'Xticklabel',[]);
    
    %%% Save figure as PNG
    % If there are multiple heatmaps in the matrix, save vascular heatmap
    if size(heatmaps,3) > 1
        fout = append(fname, '_', pathology{d});
    % Otherwise, save the pathology heatmap
    else
        fout = fname;
    end
    % If the colorbar is reversed then add suffix to filename
    if flip_cbar
        fout = append(fout, '_flip_cbar');
    end
    % Save figure as PNG
    fout = fullfile(dpath, fout);
    pause(0.1)
    saveas(gca, fout,'png');
    close;
end
end

%% Scatterplot of the pathology vs. vasculature
function lin_reg_plot(pair,sub,tstr,xl,yl,dpath,fname)
% LIN_REG_PLOT scatterplot of pathology vs. vasculature
% INPUT
%   pair (double matrix [x,2]): ROI pairs,
%                               column 1 = vasculature
%                               column 2 = pathology
%   sub (string): subject ID
%   tstr (string): figure title
%   xl (string): xlabel string
%   yl (string): ylabel string
%   dpath (string): data dicrectory path
%   fname (string): name of figure to save

% Vascular metric
x = pair(:,1);
% Pathology 
y = pair(:,2);
% Rescale Pathology from uint8 to [0,1]
y = rescale(y,0,1,'InputMin',0,'InputMax',255);
% Fit linear regression and plot
fh = figure(); fh.WindowState = 'maximized';
scatter(x,y,50,'k','filled','o');
% Axis limits, labels, etc.
ylim([0,max(y)]);
xlabel(xl);
ylabel(yl);
set(gca,'FontSize',50);
title(tstr)
legend('hide')
% Save the figure
fout = fullfile(dpath,sub,fname);
pause(0.1)
saveas(gcf,fout,'png');
close;

end