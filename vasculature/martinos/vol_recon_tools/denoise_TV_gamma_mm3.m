%% denoise_TV_gamma_mm
%% Author: Divya Varadarajan
% This function denoises the MRI using a Rician based optimization framework.
%% Inputs
% y : Noisy DWI data for one voxel (NxV), N-q-space samples(gradient directions) and V-# voxels
% lambda :  TV regularization parameter
% SNR : =1/sigma, SNR of the DWI data
% denoise_init :  flag to choose initialization. (not supported yet) -
%                 default uses noisy image as initialization

%% Outputs:
% mm_denoised_data : Denoised MRI data

%%
function [mm_denoised_data, lam] = denoise_TV_gamma_mm3(y_vol, lambda,a,b,y0,adapt_lam)

if nargin <5
    y0 =y_vol;
end

if nargin < 6
    adapt_lam=1;
end

% volume size
N1=size(y_vol,1);
N2 = size(y_vol,2);
Nslice = size(y_vol,3);
lam = lambda.*ones(N1,N2); %y/10;
options.verb=0;

for nz = 1:Nslice,
    y = y_vol(:,:,nz);
    xo = y*0;
    xnew = y0(:,:,nz);

    iter = 0;
    while (norm(xo(:)-xnew(:))/norm(xnew(:)) > 1e-3 && iter<20)
        iter = iter+1;
        xo = reshape(xnew,[N1,N2]);

        K = @(x)grad(x);
        KS = @(x)-div(x);
        Amplitude = @(u)sqrt(sum(u.^2,3));
        F = @(u)lambda*sum(sum(Amplitude(u)));
    %     G = @(x)norm(sqrt(x)./sqrt(xo) - sqrt(y)./sqrt(x),'fro')^2;
        G = @(x)norm(x - nthroot((b/a)*xo.*y.^2,3),'fro')^2;

        Normalize = @(u)u./repmat( max(Amplitude(u),1e-10), [1 1 2] );
        ProxF = @(u,tau)repmat(perform_soft_thresholding(Amplitude(u),lam.*tau), [1 1 2]).*Normalize(u);
        ProxFS = compute_dual_prox(ProxF);

        ProxG = @(x,tau)(x+tau*(nthroot((b/a)*xo.*y.^2,3)))/(1+tau);

        [xnew,~] = perform_admm(xo, K,  KS, ProxFS, ProxG,options);
        xnew(xnew<0) = 0;
        xnew(y==0)=0;
        
        if adapt_lam == 1 % spatially varying lambda that scales with intensity
            sc = xnew./max(xnew(:));
            lam = 2*sc.*lambda;
        end
    end
    xnew(isnan(xnew)) = 0;
    mm_denoised_data(:,:,nz) = reshape(abs(xnew),[N1,N2]);
end

end
