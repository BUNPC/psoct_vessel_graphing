function [ agar_status ] = view_tile_aip(ParameterFile, coordinates)

load(ParameterFile)

group = ExperimentBasic.MapIndex_Tot_offset(coordinates(1),coordinates(2)) + (Parameters.SliceID-1)*ExperimentBasic.TilesPerSlice;

XPixClip        = Parameters.XPixClip;
AipGrayRange    = Parameters.AipGrayRange;
BaseFileName    = [Scan.RawDataDir filesep Scan.FilePrefix];
GraymatterThres = Parameters.Agar.gm_aip;

figure('Position',[0,0,1200, 1200]);
for i = 1:length(group)
    VolName1 = [BaseFileName sprintf('%03i', group(i)) '_aip.nii'];
    SLO_dBI = readnifti(VolName1);
    SLO_dBI = (SLO_dBI(XPixClip:end,:));

    pct = nnz(SLO_dBI>GraymatterThres)/numel(SLO_dBI);
    agar_status(1,i) = pct;
    agar_status(2,i) = std(SLO_dBI,[],'all','omitnan');
    
    subplot(3,1,1);
    imagesc(SLO_dBI); axis equal off tight; colormap(gray); caxis([AipGrayRange]);
    title(['Slice' num2str(i) ' prc=' num2str(agar_status(1,i)) ' std=' num2str(agar_status(2,i))])
    subplot(3,1,2);
    plot(agar_status(1,:)); title('prc');
    subplot(3,1,3);
    plot(agar_status(2,:)); title('std'); ylim([1 1.5]);
    pause(0.5)
end
end