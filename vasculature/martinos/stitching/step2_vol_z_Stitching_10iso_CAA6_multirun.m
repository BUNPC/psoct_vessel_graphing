%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions required: do_downsample (path : vol_recon_beta);  
% %                   imgaussian (path : frangi filter folder)

% Add path and load parameter file
addpath('/autofs/cluster/octdata2/users/Hui/tools/dg_utils/vol_recon_beta');
% ParameterFile  = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1/Parameters.mat';
% load(ParameterFile); %Parameters Mosaic3D Scan

% %%% Set mosaic parameters
% fprintf(' - Loading Experiment file...\n %s\n', Mosaic3D.Exp);
% S = whos('-file',Mosaic3D.Exp);
% idx = find(contains({S(:).name},'Experiment_Fiji'));
% load(Mosaic3D.Exp,S(idx).name);
% Experiment  = eval(S(idx).name);
% fprintf(' - %s is loaded ...\n %s\n', S(idx).name);


%% generating sliceidx for cases with multiple runs
% NOTE: we must manually verify that there are no overlapping slies. To do
% this, navigate to the directory "StitchingFiji"
% Verify that slice[final] in run001 is different from slice001 in run002.
% This can be performed visually. 
% Perform this validation for all subsequent runs. Alternatively, examine
% the script from prior steps.
%
% We must generate an array of slice indices for each run. These slice
% indices are stored in the variable "Mosaic3D.sliceidx", which is in the
% file: /autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1/Parameters.mat'
% Steps:
%   - Manually load the parameters file.
%   - Open Mosaic3D.sliceidx
%   - This is a 3xN matrix.
%       First row = slice number within run
%       second row = overall slice number across all runs
%       third row = array of ones (always)
run1 = 1:50; % 1:[ number of slices for each run ]
run2 = 1:14;
run3 = 1:14; 

%%% Processed Directory
% Note that this will vary for each subject, depending on who ran the first
% processing. Here is one example of the subfolder:
% /autofs/cluster/octdata2/users/Chao/caa/[subID]/processed/[YYYYMMDD]/[YYYYMMDD]_[run#]_processed/
% /autofs/cluster/octdata2/users/Chao/caa/[subID]/proces_run#
% look at the parameters.mat filepath structures

multirun_p.ProcDir = {'/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run1',...
    '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run2',...
    '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/process_run3'};

%%% Generate a slice indexing matri (3xN)
%    - First row = slice number within run
%    - second row = overall slice number across all runs
%    - third row = index of the variable multirun_p.ProcDir. The run number
%               may not always start at 1.
multirun_p.in = [run1 run2 run3]; % sliceid from multiple runs
multirun_p.out = 1:length(multirun_p.in);
multirun_p.run = [ones(size(run1))*1 ones(size(run2))*2 ones(size(run3))*3];
sliceidx = [multirun_p.in;multirun_p.out;multirun_p.run];

case_notes = 'need2register';

%% load registration parameters if needed
% % % to be added into the script 
% % % this section is for interrun translation
reg_run_file = '/autofs/cluster/octdata2/users/Chao/caa/caa_6/frontal/imregcorr_p.mat';

fprintf('case notes: case %s\n',case_notes)
switch case_notes
    case 'need2register'
        load(reg_run_file,'crop_idx','tran_idx')
end
%% Set Enviroment
ProcessDir = [multirun_p.ProcDir{1} filesep 'dBI'];
% ProcessDir = '/autofs/space/omega_001/users/caa/CAA26_Occipital/Process_caa26_occipital/mus';

%% generating z stitching parameters
% if cases doesnt have multiple runs, the following code that commented out
% will work just fine. If there is only one run, uncomment the next line.
% Do not run the previous section
% sliceidx = Mosaic3D.sliceidx(1,:);
modality = 'mus';
zoff        = 0; %Mosaic3D.z_parameters(1);
z_sm        = 40; %Mosaic3D.z_parameters(2);
z_s         = 0; %Mosaic3D.z_parameters(3);
z_sms       = z_sm+z_s; 
z_m         = z_sm-z_s;

ds_xyz = [1, 1, 4];

fprintf('zoff    skirt      mainbody    skirt\n')
fprintf('%0.2i      %0.2i         %0.2i          %0.2i\n',zoff, z_s,z_m,z_s)
fprintf('%0.2i + %0.2i = %0.2ipx (%ium) -> slice thickness\n',z_s,z_m,z_sm,z_sm*2.5)

%zoff    skirt      mainbody    skirt
%0       15         25          15

%skirt + mainbody = slice thickness
%15 + 25 = 40px (100um)

% load lowest number plus 1 slice (sometimes lowest number slice can have weird intensity drop)
load([ProcessDir filesep modality '_slice001.mat'],'MosaicFinal');%c

if strcmpi(modality,'mus')
    mus_aip = mean(smooth3(MosaicFinal),3);
    thresh = multithresh(mus_aip,2); % output will be [agar-gm gm-wm]
    figure;subplot(1,2,1);imshow(mus_aip,[thresh(2),prctile(mus_aip(:),98)]); title('wm');
    subplot(1,2,2);imshow(mus_aip,[thresh(1),thresh(2)]); title('gm');


    wm_mask = mus_aip>thresh(2);    mask3d = logical(MosaicFinal(1,1,:)*0 + wm_mask);
    % gm_mask = mus_aip<thresh(2) & mus_aip>thresh(1);    mask3d = logical(Id(1,1,:)*0 + gm_mask);

    tissue=mean(reshape(MosaicFinal(mask3d),[],size(MosaicFinal,3)), [1]);clear mask3d
    tissue = tissue(:).';
    figure; plot(tissue); title('tissue decay')


    % generate weights that uniform intensity level across depth
    tissue = tissue([1:z_sms]+zoff);
    s  = tissue(z_s+1)./tissue(1:z_s);        % Top Overlapping Skirt
    ms = tissue(z_s+1)./tissue(z_s+1:z_sms);  % Non-overlapping Body  &  Bottom Overlapping Skirt
    
    % blending weight at the overlap
    switch modality
        case 'dBI'; degree = 1;   
        case 'mus'; degree = 1; % 1  
    end
    w1 = s             .*linspace(0,1,z_s).^degree; w1 = w1(:).';    % Top Overlapping Skirt
    w2 = ms(z_m+1:end) .*linspace(1,0,z_s).^degree; w2 = w2(:).';    % Bottom Overlapping Skirt
    w3 = ms(1:z_m);                                 w3 = w3(:).';

elseif strcmpi(modality,'conv')
    w1 = linspace(0,1,z_s); w1 = w1(:).';    % Top Overlapping Skirt
    w2 = linspace(1,0,z_s); w2 = w2(:).';    % Bottom Overlapping Skirt
    w3 = [1:z_m]*0+1; w3 = w3(:).';
end

nb_slices   = size(sliceidx,2);
[row_pxblock, col_pxblock] = size(MosaicFinal,[1 2]);
tot_z_pix   = z_sm*nb_slices;
Ma_sz = [row_pxblock, col_pxblock, tot_z_pix];
Ma_sz = Ma_sz./ds_xyz;                                                     %<<<<<<<<<<<<<<<<<<<<<<<<

%%% template for stitched block
Ma = zeros(Ma_sz,'single');
clear MosaicFinal s ms

%% Begin stitching  
tic

% Set output directory
outdir = '/autofs/space/omega_001/users/caa/CAA26_Occipital/Process_caa26_occipital/mus_vessel';
if ~exist(outdir, 'dir')
    mkdir(outdir)
end
fprintf(' - Output directory = \n%s\n',outdir);

% Set input directory
in_dir = multirun_p.ProcDir;
for ii = 1:length(in_dir)
    fprintf(' - Input directory = \n%s\n',in_dir{ii})
end
    
%%% loop slices
for s = 1:length(sliceidx)
    s_in = sliceidx(1,s);
    s_out = sliceidx(2,s);
    s_run = sliceidx(3,s);
    fprintf('\tslice %d in run %d \tslice %d in all runs\n',s_in,s_run,s_out);
    
    indir = in_dir{s_run};

    %%% load slice
    fdata = [indir filesep modality '_slice' sprintf('%03i',s_in) '.mat'];
    load(fdata,'MosaicFinal');
    MosaicFinal = MosaicFinal(:,:,zoff+1:z_sms);

    %eval(['data' num2str(kk) '=data;']);
    ddims = size(MosaicFinal);


    switch case_notes
        case 'need2register'
        % % % to be added into the script 
        % % % this section is for interrun translation
            c1s=crop_idx(i_run,1);
            c1e=crop_idx(i_run,2);
            c2s=crop_idx(i_run,3);
            c2e=crop_idx(i_run,4);
            t1s=tran_idx(i_run,1);
            t1e=(c1e-c1s) + tran_idx(i_run,1); % crop length + translation index
            t2s=tran_idx(i_run,2);
            t2e=(c2e-c2s) + tran_idx(i_run,2);
        case 'isuniform'
            c1s=1;
            c1e=ddims(1);
            c2s=1;
            c2e=ddims(2);
            t1s=1;
            t1e=ddims(1); % crop length + translation index
            t2s=1;
            t2e=ddims(2);
    end

    %%% stitch (cases for first, last, and other slice #s)
    if s_out==1 % first_slice
        W = reshape(repmat([w1*0+1 w3 w2],          [ddims(1)*ddims(2),1]),ddims(1),ddims(2),[]);
        zr1 = 1:z_sms/ds_xyz(3);
        zr2 = 1:z_sms+zoff; %zr2 = zr1; %

        I_in = do_downsample(MosaicFinal(:,:,zr2).*W,ds_xyz);
        Ma(t1s:t1e,t2s:t2e,zr1) = I_in(c1s:c1e,c2s:c2e,:);

    elseif s_out==length(sliceidx) % last_slice
        % Add to the `t` bottom-most pixels in Ma
        W = reshape(repmat([w1 w3 ],          [ddims(1)*ddims(2),1]),ddims(1),ddims(2),[]);            
        zr1 = [1:z_sm/ds_xyz(3)] + Ma_sz(3)-z_sm/ds_xyz(3);
        zr2 = zoff+[1:z_sm]; %zr2 = [1:z_sms]; %

        I_in = do_downsample(MosaicFinal(:,:,zr2).*W,ds_xyz);
        Ma(t1s:t1e,t2s:t2e,zr1) = Ma(t1s:t1e,t2s:t2e,zr1) + I_in(c1s:c1e,c2s:c2e,:);
    else
        % For slices (first_slice+1) to (last_slice-1)
        W = reshape(repmat([w1 w3 w2],              [ddims(1)*ddims(2),1]),ddims(1),ddims(2),[]);
        zr1 = [1:z_sms/ds_xyz(3)] + (s-1)*z_sm/ds_xyz(3);
        zr2 = zoff+[1:z_sms]; %zr2 = [1:z_sms]; %

        I_in = do_downsample(MosaicFinal(:,:,zr2).*W,ds_xyz);
        Ma(t1s:t1e,t2s:t2e,zr1) = Ma(t1s:t1e,t2s:t2e,zr1) + I_in(c1s:c1e,c2s:c2e,:);
    end
    disp('Finished slice');
end

toc

%% Save outputs
%%% raw output
resolution      = [0.01 0.01 0.0025];
resolution_ds10 = resolution .* ds_xyz;
name = [outdir filesep modality '_mean_10um-iso.slice40px.nii'];
save_mri(Ma, name, resolution_ds10)
    
%%% smoothed output
name_smooth = [outdir filesep modality '_mean_10um-iso.slice40px.sigma1.nii'];
save_mri(imgaussian(Ma,1), name_smooth, resolution_ds10)

%% Save 20um iso volume
%{
disp(' - Downsampling to 20 micron iso (z)...');

sx = 1; %20um
sy = 1; %20um
sz = 8; %2.5um -> 20um
ds_xyz = [2, 2, 2]; % sz are loaded from recon.mat

[Id20] = do_downsample(Ma, ds_xyz);
resolution_ds20 = resolution_ds10 .* ds_xyz;
name = [outdir filesep modality '_mean_20um-iso.slice40px.nii'];
save_mri(Id20, name, resolution_ds20)
%}

function save_mri(I, name, res)
% Save stitched output
    disp(' - making hdr...');
    % Make Nifti and header
    colres = res(2); 
    rowres = res(1); 
    sliceres = res(3); 
    % mri.vol = I;
    mri.volres = [res(1) res(2) res(3)];
    mri.xsize = rowres;
    mri.ysize = colres;
    mri.zsize = sliceres;
    a = diag([-colres rowres sliceres 1]);
    mri.vox2ras0 = a;
    mri.volsize = size(I);
    mri.vol = I;
    % mri.vol = flip(mri.vol,1);
    MRIwrite(mri,name,'float');
    disp(' - done - ');
end