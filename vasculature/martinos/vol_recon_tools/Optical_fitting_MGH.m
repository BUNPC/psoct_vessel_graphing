function [us,ub]= Optical_fitting_MGH(I)
% Fitting code for the data shared by Hui. 
% I: volume data for 1 tile
% s_seg: slice number,for initial test, just use arbitrary number like 1
% Z_seg: tile numberfor initial test, just use arbitrary number like 1
% datapath: datapath for original data
% zf: focus depth matrix, X*Y dimension, for initial test, just use
% arbitrary number like 1

%cd(datapath);
opts = optimset('Display','off','TolFun',1e-10); % convergence tolerence
v=ones(3,3,3)./27;  % local averaging to make the results more stable
if canUseGPU(); I = gpuArray(I);end
I=convn(I,v,'same');
I=flip(I,3);  % flip the z axis
w=1.35; % sensitivity roff-off constant, w=2.2 for 5x obj, w=2.22 for 10x obj
%I=rolloff_corr(I,w);
I=10.^(I./10); % convert to linear scale
% depth to fit, tunable
fit_depth = round(size(I,3)); 
% Z step size
Z_step=2.5; %change this after asking hui
% sensitivity roll-off correction

% for flat cut, define a constant depth as start of fitting
z0=1; % same value for all tiles, all slices
% cut out the signal above start depth
d=min(fit_depth+z0-1,size(I,3));
I=I(:,:,z0:d);


    % Average attenuation for the full ROI
    mean_I = squeeze(mean(I,[1,2]));
    mean_I = mean_I - mean(mean_I(end-5:end));
%     mean_I=flip(mean_I);
    [m,x]=max(mean_I);

    ydata1=mean_I(x:end);
    [n,x0]=find(ydata1<0.01*m);  zind=n(1);
    ydata = double(mean_I(x:zind)')/m;
%     z = (x:(length(ydata)+x-1))*Z_step;
    z = (x:zind)*Z_step;
%     
    fun = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*zdata).*(1./(1+((zdata-p(3))./p(4)).^2)));
    lb = [0.1 0.005 2*Z_step 1]; %ask Hui what is the expected rayleigh range
    ub = [100 0.02 80*Z_step 200];
    
    if canUseGPU(); [z,ydata] = gather(z,ydata);end
    est = lsqcurvefit(fun,[4 0.005 20*Z_step 80],z,ydata,lb,ub,opts);
%     
     A = fun(est,z);
      %plotting intial fitting - commented out when running large batches
        figure
        plot(z,ydata,'b.')
        hold on
        plot(z,A,'r-')
        xlabel('z (um)')
        ylabel('I')
        title('Four parameter fit of averaged data')
        dim = [0.2 0.2 0.3 0.3];
        str = {'Estimated values: ',['Relative back scattering: ',num2str(est(1),4)],['Scattering coefficient: ',...
            num2str(est(2)*1000,4),'mm^-^1'],['Focus depth: ',num2str(est(3),4),'um'],['Rayleigh estimate: ',num2str(round(est(4)),4),'um']};
        annotation('textbox',dim,'String',str,'FitBoxToText','on');

    %% Curve fitting for the whole image
    res = 1;  %A-line averaging factor, 5 means 10x10 averaging
    est_pix = zeros(round(size(I,2)/res-1),round(size(I,3)/res-1),3);
%         est_pix = zeros(round(size(I,2)/res-1),round(size(I,3)/res-1),2);

    for i = 1:round(size(I,1)/res)
        for j = 1:round(size(I,2)/res)
            imin=min((i+1)*res,size(I,1));
            jmin=min((j+1)*res,size(I,2));
            area = I((i-1)*res+1:imin,(j-1)*res+1:jmin,:);
            int = squeeze(mean(area,[1,2]));
            int = (int - mean(int(end-5:end)));
%             xloc=findchangepts(int);
            m=max(int);
            xloc=x;   % empirical start point of fitting is good enough for flat slice
%             [n,x0]=find(int<0.05*m);  zind=x0(1);
    
            if m > 0 % change this threshold to avoid fitting agar and save some time
                l=min(size(int,1),xloc+fit_depth-1);
                z = (xloc:(length(ydata)+xloc-1))*Z_step;
                ydata_tile = double(int(xloc:(length(ydata)+xloc-1))')/m;

%                 l=zind;
%                 ydata = double(int(xloc:l)');
%                 z = (xloc:l)*Z_step;

%                 fun = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*zdata).*(1./(1+((zdata-p(3))./p(4)).^2)));
%                 lb = [1e8 0.001 2*Z_step 10]; %ask Hui what is the expected rayleigh range
%                 ub = [1e20 0.05 80*Z_step 500];
%                 est_pix(i,j,:) = lsqcurvefit(fun,[1e10 0.01 20*Z_step 80],z,ydata_tile,lb,ub,opts);
    
                fun_pix = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*(zdata)).*(1./(1+((zdata-p(3))./(est(4))).^2))); % 3-parameter fitting using empirical rayleigh range
                lb = [0.1 0.0005 2*Z_step 1]; 
                ub = [100 0.02 80*Z_step 200];
%                 temp = lsqcurvefit(fun_pix,[4 0.005 20*Z_step],z,ydata_tile,lb,ub,opts);
                if canUseGPU(); [z,ydata_tile] = gather(z,ydata_tile);end
                est_pix(i,j,:) = lsqcurvefit(fun_pix,[4 0.005 20*Z_step],z,ydata_tile,lb,ub,opts);
                
%                 fun_pix = @(p,zdata)sqrt(p(1).*exp(-2.*p(2).*(zdata)).*(1./(1+((zdata-est(3))./(est(4))).^2))); % 2-parameter fitting using tile-dependent zf and rayleigh range from above fitting
%                 % upper and lower bound for fitting parameters
%                 lb = [0.1 0.0005];     ub=[100 0.02];
%                 temp= lsqcurvefit(fun_pix,[4 0.005],z,ydata_tile,lb,ub,opts);
% %                 est_pix(i,j,:) = lsqcurvefit(fun_pix,[4 0.005],z,ydata_tile,lb,ub,opts);
%                 
%                 A = fun_pix(temp,z);
%                 figure
%                 plot(z,ydata_tile,'b.')
%                 hold on
%                 plot(z,A,'r-')
                
            else
                est_pix(i,j,:) = [0 0 0];
            end
            
        end
    end

    us = 1000.*squeeze(est_pix(:,:,2));     % unit:mm-1
    ub = squeeze(est_pix(:,:,1));
