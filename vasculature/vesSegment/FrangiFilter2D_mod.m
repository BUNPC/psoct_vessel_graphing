function outIm = FrangiFilter2D(I, sigma, ci, bg)
% This function FRANGIFILTER2D uses the eigenvectors of the Hessian to
% compute the likeliness of an image region to vessels, according
% to the method described by Frangi:2001 (Chapter 2).
%
% J = FrangiFilter2D(I, Options)
%
% inputs,
%   I : The input image (vessel image)
%   sigma : The hessian filter size (two-element vector for upper and lower
%   bounds)
%   ci : scaling parameter for vesselness filter, default : ci = 15;
%   bg : bright background (1) or dark backgraound (0), default : bg =1
%
% outputs,
%   J : The vessel enhanced image (pixel is the maximum found in all scales) 
%
% Modified by Jiarui Yang

    % make sure we have all input arguments
    if nargin < 2
        sigma=1;
        ci=15;
        bg=1;
    elseif nargin < 3
        ci=15;
        bg=1;
    elseif nargin < 4
        bg=1;
    end
    
    sigmas=sigma(1):sigma(2);
    sigmas = sort(sigmas, 'ascend');

    beta  = 2*(0.5^2);
    c     = 2*(ci^2);

    % Make matrices to store all filterd images
    ALLfiltered=zeros([size(I) length(sigmas)]);
    % ALLangles=zeros([size(I) length(sigmas)]);

    h = waitbar(0,'Please wait... performing 2D fiber enhancement');

    % Frangi filter for all sigmas
    for i = 1:length(sigmas)
        % Show progress
        waitbar((i-1)/length(sigma));

        % Make 2D hessian
        [Dxx,Dxy,Dyy] = Hessian2D(I,sigmas(i));

        % Correct for scale
        Dxx = (sigmas(i)^2)*Dxx;
        Dxy = (sigmas(i)^2)*Dxy;
        Dyy = (sigmas(i)^2)*Dyy;

        % Calculate (abs sorted) eigenvalues and vectors
        [Lambda2,Lambda1,~,~]=eig2image(Dxx,Dxy,Dyy);

        % Compute the direction of the minor eigenvector
        % angles = atan2(Ix,Iy);

        % Compute some similarity measures
        Lambda1(Lambda1==0) = eps;
        Rb = (Lambda2./Lambda1).^2;
        S2 = Lambda1.^2 + Lambda2.^2;

        % Compute the output image
        Ifiltered = exp(-Rb/beta) .*(ones(size(I))-exp(-S2/c));

        % see pp. 45
        if(bg==1)
            Ifiltered(Lambda1<0)=0;
        else
            Ifiltered(Lambda1>0)=0;
        end
        % store the results in 3D matrices
        ALLfiltered(:,:,i) = Ifiltered;
        % ALLangles(:,:,i) = angles;
    end
    close(h);
    
    % Return for every pixel the value of the scale(sigma) with the maximum 
    % output pixel value
    if length(sigmas) > 1
        [outIm, ~] = max(ALLfiltered,[],3);
        outIm = reshape(outIm,size(I));
    else
        outIm = reshape(ALLfiltered,size(I));
    end
end