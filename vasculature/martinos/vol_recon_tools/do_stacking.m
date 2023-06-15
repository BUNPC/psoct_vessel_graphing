function [] = do_stacking(ParameterFile, modality)

load(ParameterFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING SLICE INPUT 
sliceid     = Stack.sliceid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING DIRECTORIES 
indir       = Stack.indir;
outdir      = Stack.outdir; if ~exist(outdir,'dir'); mkdir(outdir); end

fprintf(' - Input directory = \n%s\n',indir);
fprintf(' - Output directory = \n%s\n',outdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SETTING RESOLUTION 
resolution = Stack.Resolution;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Begin!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Set Waitbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = waitbar(0,'1','Name',['Stacking ' modality ' ...'],...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

%% Stacking
for i = 1:length(sliceid)
    
    Imag = [indir filesep  modality '_slice' sprintf('%03i',sliceid(i)) '.mat'];
    load(Imag);
    S = whos('-file',Imag);
    I = eval(S.name);
    if ~isa(I,'single');single(I);end
    
    if i == 1; stacknii = zeros([size(I,[1,2]),length(sliceid)],'like',I) ;end %fprintf('Stacking %s ',modality)
    
    stacknii(:,:,i) = I;
    
    % Update waitbar and message
    step = i; steps = length(sliceid);
    % Check for clicked Cancel button
    if getappdata(f,'canceling')
        break
    end
    waitbar(step/steps,f,sprintf('Slice %i/ -%i- /%i',sliceid(1), sliceid(i), sliceid(end)))
end
stacknii = rot90(stacknii);
%% Saving

waitbar(step/steps,f,sprintf('Saving Stacked %s.nii', modality))

name = [outdir filesep sprintf('Stacked_%s.nii', modality)];

fprintf(' xy res=%g mm\n z  res=%g mm\n',resolution(1),resolution(3));


disp(' - making hdr...');
% Make Nifti and header
colres = resolution(2); 
rowres = resolution(1); 
sliceres = resolution(3); 
% mri.vol = I;
mri.volres = [resolution(1) resolution(2) resolution(3)];
mri.xsize = rowres;
mri.ysize = colres;
mri.zsize = sliceres;
a = diag([-colres rowres sliceres 1]);

mri.vox2ras0 = a;
mri.volsize = size(stacknii);
mri.vol = stacknii;
% mri.vol = flip(mri.vol,1);
MRIwrite(mri,name,'float');
disp(' - END - ');

delete(f)