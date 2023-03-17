%% Curve fitting
clear all;

datapath  = strcat('C:\Users\jryang\Desktop\Data\191029phantom\20\depth14\');

% add subfunctions for the script
addpath('C:\Users\jryang\Desktop\Data\code\');

tic
% get the directory of all image tiles
cd(datapath);
%% load data
load(strcat('avg.mat'));

opts = optimset('Display','off','TolFun',1e-10);        % fitting options
%% FOV curvature correction, sensitivity roll-off correction and volumetric averaging

    %s = data_ori;39.8
%     n=3;
%     v=ones(3,n,n)./(n*n*3);
    %load(filename0(iFile).name);
   %    s = curvature_corr(s);    % field curvature correction
%    s = convn(s, v, 'same');
    %s = rolloff_corr(s, 2.3);      % sensitivity rolloff correction

%% Fit for averaged A-line
    avg = squeeze(mean(mean(s, 2), 3)); % calculate averaged aline
    avg = avg - mean(avg(end-19:end));                   % remove noise floor
%     if depth < 3
%         start_index = 21; 
%     else
%         start_index = depth*25;
%     end
    fit_depth = round(390/3);       % depth want to fit, tunable
    start_index = 71;
    [m, x] = max(avg(start_index:end));                           % depth of first pixel in fitting range
    x = x + start_index - 10;
    ydata = avg(x:x+fit_depth-1);
    ydata = double(ydata');
    z = x*3:3:(length(ydata)*3)+((x-1)*3);
    fun = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*zdata).*(1./(1+((zdata-p(3))./p(4)).^2)));
    lb = [0 0 (x-20)*3 0];
    ub = [30000 0.01 (x+200)*3 80];
    [est, residual] = lsqcurvefit(fun,[100 0.0005 (x+50)*3 70],z,ydata,lb,ub,opts);
%    A = fun(est,z);
    % plotting intial fitting
%     figure;
%     plot(z,ydata,'b.');
%     hold on;
%     plot(z,A,'r-');
%     xlabel('z (um)');
%     ylabel('I');
%     title('Four parameter fit of averaged data');
%     dim = [0.2 0.2 0.3 0.3];
%     str = {'Estimated values: ',['Relative back scattering: ',num2str(est(1),4)],['Scattering coefficient: ',...
%         num2str(est(2)*1000,4),'mm^-^1'],['Focus depth: ',num2str(est(3),4),'um'],['Rayleigh estimate: ',num2str(round(est(4)),4),'um']};
%     annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
%% Curve fitting for the whole image, using effective rayleigh range
res = 20;
est_pix = zeros(round(size(s,2)/res),round(size(s,3)/res),4);

for i = 1:round(size(s,2)/res)
    for j = 1:round(size(s,3)/res)
        temp = squeeze(mean(mean(s(:,(i-1)*res+1:i*res,(j-1)*res+1:j*res),2),3));
        temp = temp - mean(temp(end-19:end));
        [m,x] = max(temp(start_index:end));
        x = x + start_index - 10;
        ydata = temp(x:x+fit_depth-1);
        ydata = double(ydata');
        z = x*3:3:(length(ydata)*3)+((x-1)*3);
        fun_pix = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*zdata).*(1./(1+((zdata-p(3))./p(4)).^2)));
        lb = [0 0 est(3)-100 0];
        ub=[30000 0.01 est(3)+100 80];
        est2 = lsqcurvefit(fun_pix,[100 0.001 est(3) est(4)],z,ydata,lb,ub,opts);
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
%             num2str(est2(2)*1000,4),'mm^-^1'],['Focus depth: ',num2str(est2(3),4),'um'],['Rayleigh estimate: ',num2str(round(est(4)),4),'um']};
%         annotation('textbox',dim,'String',str,'FitBoxToText','on');
%         disp(strcat('B-scan No.',num2str(i),' is finished'));
%     end
end
%info = strsplit(filename,'-');
%z_seg = str2double(info{2}(1:3));
% z_seg(2) = str2double(info{2}(7:9));
%savename='ScatFit';
%save([savename, '.mat'],'est_pix')



%% visualization

%us = 1000.*squeeze(est_pix(:,:,2));
% figure;
% imagesc((1:20)*0.05,(1:20)*0.05,us);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% xticklabels({0 0.2 0.4 0.6 0.8 1});
% yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({0 0.2 0.4 0.6 0.8 1});
% set(gca,'TickLength',[0 0]);
% caxis([prctile(us(:),5),prctile(us(:),95)]);
% colorbar;
% xlabel('x (mm)');
% ylabel('y (mm)');
% title('Scattering coefficient (mm^-^1)');
% 
% ub = squeeze(est_pix(:,:,1));
% figure;
% imagesc((1:50)*0.02,(1:50)*0.02,ub);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% xticklabels({0 0.2 0.4 0.6 0.8 1});
% yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({0 0.2 0.4 0.6 0.8 1});
% set(gca,'TickLength',[0 0]);
% caxis([prctile(ub(:),5),prctile(ub(:),95)]);
% xlabel('x (mm)');
% ylabel('y (mm)');
% title('Relative back scatter');
% colorbar;
% 
% zf = squeeze(est_pix(:,:,3));
% zf = zf - min(zf(:));
% h1=figure;
% imagesc((1:20)*0.05,(1:20)*0.05,zf);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% xticklabels({0 0.2 0.4 0.6 0.8 1});
% yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({0 0.2 0.4 0.6 0.8 1});
% set(gca,'TickLength',[0 0]);
% caxis([9 34]);%([floor(prctile(zf(:),5)),floor(prctile(zf(:),95))]);
% xlabel('x (mm)');
% ylabel('y (mm)');
% %title('depth of focus (um)');
% cb1=colorbar('h1');
% set(cb1,'Location','eastoutside');
% set(cb1,'YTick',9:5:34);%floor(prctile(zf(:),5)):5:floor(prctile(zf(:),95)));

% h1=figure;
% imagesc((1:50)*0.02,(1:50)*0.02,sur);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% xticklabels({0 0.2 0.4 0.6 0.8 1});
% yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({0 0.2 0.4 0.6 0.8 1});
% set(gca,'TickLength',[0 0]);
% caxis([floor(prctile(sur(:),5)),floor(prctile(sur(:),95))]);
% xlabel('x (mm)');
% ylabel('y (mm)');
% %title('depth of focus (um)');
% cb1=colorbar('h1');
% set(cb1,'Location','eastoutside');
% set(cb1,'YTick',floor(prctile(sur(:),5)):floor(prctile(sur(:),95)));

zrs = squeeze(est_pix(:,:,4));
% mean2(zrs)
% std2(zrs)
 h2=figure;
% imagesc((1:50)*0.02,(1:50)*0.02,zrs);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% xticklabels({0 0.2 0.4 0.6 0.8 1});
% yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({0 0.2 0.4 0.6 0.8 1});
% set(gca,'TickLength',[0 0]);
% caxis([floor(prctile(zrs(:),5)) floor(prctile(zrs(:),95))]);
% xlabel('x (mm)');
% ylabel('y (mm)');
% %title('Rayleigh range (um)');
% cb2=colorbar('h2');
% set(cb2,'Location','eastoutside');
% set(cb2,'YTick',floor(prctile(zrs(:),5)):floor(prctile(zrs(:),95)));
imagesc((1:20)*0.05,(1:20)*0.05,zrs);
xticks([0 0.2 0.4 0.6 0.8 1]);
xticklabels({0 0.2 0.4 0.6 0.8 1});
yticks([0 0.2 0.4 0.6 0.8 1]);
yticklabels({0 0.2 0.4 0.6 0.8 1});
set(gca,'TickLength',[0 0]);
caxis([floor(prctile(zrs(:),5)) floor(prctile(zrs(:),95))]);
xlabel('x (mm)');
ylabel('y (mm)');
%title('Rayleigh range (um)');
cb2=colorbar('h2');
set(cb2,'Location','eastoutside');
set(cb2,'YTick',floor(prctile(zrs(:),5)):floor(prctile(zrs(:),95)));
% 
% aip = squeeze(mean(s(x:x+fit_depth-1,:,:),1));
% figure;
% imagesc((1:400)*0.0025,(1:400)*0.0025,aip);
% xticks([0 0.2 0.4 0.6 0.8 1]);
% xticklabels({0 0.2 0.4 0.6 0.8 1});
% yticks([0 0.2 0.4 0.6 0.8 1]);
% yticklabels({0 0.2 0.4 0.6 0.8 1});
% set(gca,'TickLength',[0 0]);
% caxis([prctile(aip(:),5),prctile(aip(:),95)]);
% xlabel('x (mm)');
% ylabel('y (mm)');
% title('AIP');
% colorbar;

toc
