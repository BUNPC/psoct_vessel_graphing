clear all;
FileNum = 2:10;
zf_true0 = (0:50:490);
zf_true = zf_true0(FileNum-1);

load('/space/vault/3/users/Hui/ThorlabsCalibration_101416/sensitivity_10x/rolloff');
Rolloff = (rolloff); % 1.1; 1.1; 1.1; 0.95; 0.9; 

SubPath = 'Dil5/';
FileMus = 2;

zoff1 = 30;
for ii = 1:length(FileNum)
    load(strcat([SubPath 'dBI_', sprintf('%03i', FileNum(ii))]));
    Itemp = mean(IJones1,2)./Rolloff;
%     figure,plot(Itemp)
    [maxI,maxInd] = max(Itemp(zoff1+1:zoff1+100,:));
    zoffAll(ii) = maxInd + zoff1;
    
    IoffAll(ii) = mean(Itemp(300:320,1));
    I(:,:,ii) = IJones1./repmat(Rolloff,1,size(IJones1,2));
end
if SubPath == 'Dil4/'
    zoffAll = [109;107;108;108;110;112;108;109;108];
end


%%
figure(6)
subplot(1,2,1)
imagesc(10*log10(mean(I(:,:,1),3)));colormap gray

subplot(1,2,2)
imagesc(10*log10(mean(I(:,:,3),3)));colormap gray

% figure(2)
% plot(squeeze(mean(I(:,130:450,:),2)))


%%

n_avg = 200;
xoff = 0;
for kk = 1:floor(size(I,2)/n_avg)
    kk
    II2 = mean(I(1:end,xoff+1+(kk-1)*n_avg:xoff+(kk)*n_avg,:),2);
     
    %%
    % fit
    global g
    xfitAll = [];
    
%     figure(1)
    for iD = 1:size(I,3)
        %%
        zoff = zoffAll(iD);
        zoff_start = 15;
        % data 102116: lst = zoff+11:zoff+190;
        lst = zoff+zoff_start :zoff+zoff_start+204;
        II3(:,iD) = II2(lst,iD) - IoffAll(iD);
        x = ([1:size(II2,1)]) * 3.03; % 2.9 um per pixel
        x2 = x(lst);
        xx = ([-50:size(II2,1)]) * 3.03;
        
        
        g.zz = zf_true(iD);
        g.zs = zoff_start;
        % 'levenberg-marquardt'
        OPTIONS = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','TolFun',1e-10,'Display','none');
        
        % according to Yun... de Boer 2003
        % Ioff = (sin(eta)/eta)^2 exp(-w^2/2ln2 eta^2) where eta = pi/2 Z/Zrd
        % Zrd = lambda^2 / (4 \Delta lambda) where \Delta lambda is spacing
        % between pixels
        % w = \delta lambda / \Delta lambda where \delta lambda is spectrometer
        % resolution (i.e. FWHM).
        lx = size(II3,1);
        if 0
            xoffIdx = [Surface(nx,ny)+[-20:-10 (zmax-10):zmax]]';
            xoff = (xoffIdx-Surface(nx,ny)) * 2.9;
            yoff = mean(mean( I(xoffIdx,nx+[0:(havg-1)], ny+[0:(havg-1)]),2),3);
            p = polyfit(xoff,yoff,1);
            Ioff = p(2);
            IoffS = p(1);
        elseif 0
            Ioff = 100;
            IoffS = 0.01;
        end
        
        nParam =3;
        x0 = [4 0 78];%David attenuation (1/mm), focus depth (um), D (um)
        
        yy = II3(:,iD)/1e6;
        x = ([1:size(II3,1)]) * 3.03; % 2.9 um per pixel
        
        g.yy = yy;
        %    g.Ioff = Ioff;
        %    g.IoffS = IoffS;
        
        [xfit,resnorm,res] = lsqcurvefit(@function_Profile2, x0, x, yy',[-inf -inf -inf],[inf inf inf],OPTIONS );
%         [xfit,resnorm,res] = lsqcurvefit(@function_Profile2, x0, x, yy',[-5 -200 -100],[20 200 300],OPTIONS );
        xfitAll(1:nParam,iD) = xfit(1:nParam);
        xfitAll(nParam+1,iD) = g.Io;
        xfitAll(nParam+2,iD) = resnorm;
        xfitAll(nParam+3,iD) = 1-resnorm/sum((yy-mean(yy)).^2);
        
        %yyOff = Ioff+IoffS*xoff;%(xAll-Surface(nx,ny));
        
        if mod(iD,2)==0
        subplot(1,4,ceil(iD/2))
        %plot(x,yy,'b.',xoff,(yyOff),'b-',x,(yy),'ko',x,(function_Profile(xfit(1:nParam),x)),'k-')
        plot(x,yy,'b.',x,(function_Profile2(xfit(1:nParam),x)),'k-','linewidth',1.2)
        set(gca,'fontsize',16,'xtick',0:200:600,'xticklabel',35:200:635);
%         title( sprintf('zf=%d \mum', round(40*1.33*((iD-1)/2+2))) )
        title( ['zf=', num2str(round(50*(iD-1))), '\mum'] )
        ylim([0 7.5])
        xlabel('depth (\mum)','fontsize',16);
        ylabel('I (A.U.)','fontsize',16);
        xlim([0 600])
        end
        
%         subplot(3,5,iD)
%         %plot(x,yy,'b.',xoff,(yyOff),'b-',x,(yy),'ko',x,(function_Profile(xfit(1:nParam),x)),'k-')
%         plot(xx,(function_Profile2(xfit(1:nParam),xx)),'k-')
    end
%     figure(11)
%     subplot(1,2,1)
%     semilogy(II3)
%     
%     subplot(1,2,2)
%     semilogy(x2,II3)
   
    xfitAll2(:,:,kk) = xfitAll;
end


% %%
figure(10)
y_mus = mean(xfitAll2(1,:,:),3);
err_mus = std(xfitAll2(1,:,:),[],3);%/sqrt(size(xfitAll2,3));
errorbar(zf_true-zf_true(1),y_mus,err_mus,'o-')
set(gca,'fontsize',14)
xlim([-5 440]);
ylim([-5 20])
title('Obj: 10x, Dil: 5x, mus')
xlabel('focus depth')
ylabel('mus (1/mm)')
% sc=err_mus./y_mus

figure(20)
y_zr = median(xfitAll2(3,:,:),3);
err_zr = std(xfitAll2(3,:,:),[],3)/sqrt(size(xfitAll2,3));
errorbar(zf_true-zf_true(1),y_zr,err_zr,'o-')
set(gca,'fontsize',14)
xlim([-5 430]);
% ylim([0 350])
title('Obj: 5x, Dil: 5x, Rayleigh range')
xlabel('focus depth')
ylabel('zR (um)')

figure(30)
y_zf = median(xfitAll2(2,:,:),3);
err_zf = std(xfitAll2(2,:,:),[],3)/sqrt(size(xfitAll2,3));
errorbar(zf_true-zf_true(1),y_zf,err_zf,'o-')
set(gca,'fontsize',14)
xlim([-5 430]);
% ylim([-300 600])
title('Obj: 5x, Dil: 5x, estimated focus depth')
xlabel('focus depth')

figure(40)
y_I0 = median(xfitAll2(4,:,:),3);
err_I0 = std(xfitAll2(4,:,:),[],3)/sqrt(size(xfitAll2,3));
errorbar(zf_true-zf_true(1),y_I0,err_I0,'o-')
set(gca,'fontsize',14)
xlim([-5 430]);
% ylim([-300 600])
title('Obj: 5x, Dil: 5x, I0')
xlabel('I0')
% % 
% % 
% musCon = squeeze(xfitAll2(1,:,:));
% save(strcat([SubPath 'mus_20170525']),'musCon');

zfa = squeeze(xfitAll2(2,:,:));
save(strcat([SubPath 'zf_20170525']),'zfa');
zrsa = squeeze(xfitAll2(3,:,:));
save(strcat([SubPath 'zrs_20170525']),'zrsa');

% 
% Cparam4=[];
% for ii = 1:length(FileNum)
% paramatrix = (squeeze(xfitAll2(1:4,ii,:)))';
% Cparam4 = [Cparam4; cov(paramatrix)];
% end
% 
% 
% figure,plot([20,40,80,160,320,640],sc,'*','markersize',8,'linewidth',2)
% hold on,plot(20:640,1./sqrt(20:640),'color',[0 0 0],'linewidth',2)
% set(gca,'fontsize',16)
% legend 'zf1' 'zf2' 'zf3' 'zf4' 'zf5' 'zf6' 'zf7' 'zf8' 'zf9' 'theoretical'
% xlabel('# of average','fontsize',16);
% ylabel('speckle contrast','fontsize',16);
% ylim([0 1])
% xlim([0 650])
% title('mus (4/mm)','fontsize',16);
% 
% % sc2=sc./repmat([25.17;16.09;6.978;1.767;1.147;0.9072;0.7562;0.6652;0.8188],10,1);
% figure,plot([1,2,5,10,20,40,80,160,320,640],sc2,'o-','linewidth',1)
% hold on,plot(1:640,1./sqrt(1:640),'color',[0 0 0],'linewidth',2)
% set(gca,'fontsize',16)
% legend 'zf1' 'zf2' 'zf3' 'zf4' 'zf5' 'zf6' 'zf7' 'zf8' 'zf9' 'theoretical'
% xlabel('# of average','fontsize',16);
% ylabel('normalized StD','fontsize',16);
% ylim([0 1])
% xlim([0 650])
% title('mus (8/mm)','fontsize',16);



% for ii = 4%1:length(FileNum)
% paramatrix(:,1) = (squeeze(xfitAll2(1,ii,:)))';
% paramatrix(:,2) = (squeeze(xfitAll2(4,ii,:)))';
% paramatrix(:,3) = (squeeze(xfitAll2(2,ii,:)))';
% paramatrix(:,4) = (squeeze(xfitAll2(3,ii,:)))';
% save(['CovAnalysis/20170525/param4_mus' num2str(FileMus) '_zf' num2str(FileNum(ii))]);
% end

