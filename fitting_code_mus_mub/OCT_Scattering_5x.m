%% Curve fitting - estimate effective Reylaigh range and focus depth
clear all;
datapath  = strcat('C:\Users\jryang\Desktop\30-s3\');

% add subfunctions for the script
addpath('C:\Users\jryang\Desktop\Data\code\');

opts = optimset('Display','off','TolFun',1e-10);        % fitting options

% get the directory of all image tiles
cd(datapath);
load('data.mat');  

%     n=3;
%     v=ones(3,n,n)./(n*n*3);
    %load(filename0(iFile).name);
  %  s = curvature_corr(s);    % field curvature correction
%    s = convn(s, v, 'same');
    %s = rolloff_corr(s, 2.3);      % sensitivity rolloff correction

    avg = squeeze(mean(mean(s, 2), 3)); % calculate averaged aline
    avg = avg - mean(avg(end-19:end));                   % remove noise floor

    fit_depth = round(450/3);       % depth want to fit, tunable
    [m, xav] = max(avg(31:end));                           % depth of first pixel in fitting range
    xav = xav + 130;
    ydata = avg(xav:xav+fit_depth-1);
    ydata = double(ydata');
    z = xav*3:3:(length(ydata)*3)+((xav-1)*3);
    fun = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*zdata));
    lb = [0 0];
    ub = [100 0.01];
    est = lsqcurvefit(fun,[10 0.001],z,ydata,lb,ub,opts);
   A = fun(est,z);
   % plotting intial fitting
    figure;
    plot(z,ydata,'b.');
    hold on;
    plot(z,A,'r-');
    xlabel('z (um)');
    ylabel('I');
    title('Four parameter fit of averaged data');
    dim = [0.2 0.2 0.3 0.3];
    str = {'Estimated values: ',['Relative back scattering: ',num2str(est(1),4)],['Scattering coefficient: ',...
        num2str(est(2)*1000,4),'mm^-^1']};
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
%% Curve fitting for the whole image, using effective rayleigh range
%tic
res = 10;
est_pix = zeros(round(size(s,2)/res),round(size(s,3)/res),2);

for i = 1:round(size(s,2)/res)
    for j = 1:round(size(s,3)/res)
        temp = squeeze(mean(mean(s(:,(i-1)*res+1:i*res,(j-1)*res+1:j*res),2),3));
        temp = movmean(temp,5);
        temp = temp - mean(temp(350:end));
        [m,x] = max(temp(21:200));
        x = x + 30;
        ydata = temp(x:x+fit_depth-1);
        ydata = ydata';
        z = x*3:3:(length(ydata)*3)+((x-1)*3);
        fun_pix = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*zdata));
        lb = [0 0];
        ub=[10000 0.05];
        est2 = lsqcurvefit(fun_pix,[100 0.001],z,ydata,lb,ub,opts);
        est_pix(i,j,:) = est2;
    end
%     if(i==15)
%         A = fun_pix(est2 ,z);
%         % plotting intial fitting
%         figure;
%         plot(z,ydata,'b.');
%         hold on;
%         plot(z,A,'r-');
%         xlabel('z (um)');
%         ylabel('I');
%         title('Four parameter fit of averaged data');
%         dim = [0.2 0.2 0.3 0.3];
%         str = {'Estimated values: ',['Relative back scattering: ',num2str(est2(1),4)],['Scattering coefficient: ',...
%             num2str(est2(2)*1000,4),'mm^-^1']};
%         annotation('textbox',dim,'String',str,'FitBoxToText','on');       
%     end
%    disp(strcat('B-scan No.',num2str(i),' is finished'));
end
%info = strsplit(filename,'-');
%z_seg = str2double(info{2}(1:3));
% z_seg(2) = str2double(info{2}(7:9));
% savename='ScatFit';
% save([savename, '.mat'],'est_pix')
us=1000.*squeeze(est_pix(:,:,2));


%% visualization

% us = 1000.*squeeze(est_pix(:,:,2));
% figure;
% imagesc((1:50)*0.02,(1:50)*0.02,us)
% caxis([prctile(us(:),5),prctile(us(:),95)]);
% colorbar
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('Scattering coefficient (mm^-^1)')
% 
% ub = squeeze(est_pix(:,:,1));
% figure;
% imagesc((1:50)*0.02,(1:50)*0.02,ub)
% caxis([prctile(ub(:),5),prctile(ub(:),95)]);
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('Relative back scatter')
% colorbar
% 
% zf = squeeze(est_pix(:,:,3));
% figure;
% imagesc((1:50)*0.02,(1:50)*0.02,zf)
% caxis([prctile(zf(:),5),prctile(zf(:),95)]);
% %caxis([240,300]);
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('depth of focus (um)')
% colorbar

% zrs = squeeze(est_pix(:,:,4));
% figure;
% imagesc((1:50)*0.02,(1:50)*0.02,zrs)
% caxis([prctile(zrs(:),5),prctile(zrs(:),95)]);
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('Rayleigh range (um)')
% colorbar

% aip = squeeze(mean(s(xav:xav+fit_depth-1,:,:),1));
% figure;
% imagesc((1:400)*0.0025,(1:400)*0.0025,aip)
% caxis([prctile(aip(:),5),prctile(aip(:),95)]);
% xlabel('x (mm)')
% ylabel('y (mm)')
% title('AIP')
% colorbar
% 
% toc