%% Curve fitting - estimate effective Reylaigh range and focus depth
clear all;
datapath  = strcat('C:\Users\jryang\Desktop\Data\0629TDE\sca_length\day0\');

% add subfunctions for the script
addpath('C:\Users\jryang\Desktop\Data\code\');

% get the directory of all image tiles
cd(datapath);
filename0=dir('*.dat');
for iFile=1:length(filename0)    
    %% get data information
    nk = 400; nxRpt = 1; nx=400; nyRpt = 1; ny = 400;
    dim=[nk nxRpt nx nyRpt ny];

    %% image reconstruction
    ifilePath=[datapath,filename0(iFile).name];
    disp(['Start loading file ', datestr(now,'DD:HH:MM')]);
    [data_ori] = ReadDat_single(ifilePath, dim); 
% 
    s = data_ori;
%     n=3;
%     v=ones(3,n,n)./(n*n*3);
    %load(filename0(iFile).name);
    s = curvature_corr(s);    % field curvature correction
%    s = convn(s, v, 'same');
    %s = rolloff_corr(s, 2.3);      % sensitivity rolloff correction

    avg = squeeze(mean(mean(s, 2), 3)); % calculate averaged aline
    avg = avg - mean(avg(350:end));                   % remove noise floor
    avg (avg(:)<=0) = 0.00001;
    avg = log(movmean(avg,5));

    fit_depth = round(30/3);       % depth want to fit, tunable
    [m, xav] = max(avg(21:end));                           % depth of first pixel in fitting range
    xav = xav + 20;
    xav = 40;
    ydata = avg(xav:xav+fit_depth-1);
    ydata = ydata';
    z = xav*3:3:(length(ydata)*3)+((xav-1)*3);
    est = polyfit(z,ydata,1);
    A = est(1)*z+est(2);
    mus = (-est(1)/2)*1000;        % in mm^-1
    mub = exp(est(2));
    % plotting intial fitting
    figure;
    plot(z,ydata,'b.');
    hold on;
    plot(z,A,'r-');
    xlabel('z (um)');
    ylabel('I');
    title('Linear fit of scattering length');
    dim = [0.2 0.2 0.3 0.3];
    str = {'Estimated values: ',['slope: ',num2str(est(1),4)],['intercept: ',...
        num2str(est(2),4)],['Scattering coefficient: ',num2str((mus),4),'mm^-^1'],...
        ['Backscattering coefficient: ',num2str((mub),4)]};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    


end