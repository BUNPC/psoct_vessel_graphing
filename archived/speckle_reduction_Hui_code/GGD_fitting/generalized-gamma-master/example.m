img=squeeze(slice(55,1:500,501:1000));
img=(img-min(img(:)))./(max(img(:))-min(img(:)));
img(img(:)<=0)=0.001;
img=vect(img);
%thresh=prctile(img,5);
%img(img(:)<=thresh)=[];
img=img./mean(img(:));
pdEstimated = GeneralizedGamma();
param = pdEstimated.fitDist(img);
% compare with theoretical model
pdTrue = GeneralizedGamma(param(1),param(2),param(3));
sample = pdTrue.drawSample(1000000);
% plot
edges = 0:0.01:3;
h1 = histcounts(img,edges);
h2 = histcounts(sample,edges);
h1=h1./sum(h1);
h2=h2./sum(h2);

f=figure;
% set(f,'position',[100 100 1200 400]);
a1=area(edges(2:end),h1);
xlabel('Normalized Intensity (au)');
ylabel('Fraction of total samples');
a1.FaceAlpha=0.5;
hold on;
a2=area(edges(2:end),h2);
a2.FaceAlpha=0.5;
legend('Actual data','Drawn samples');
set(gca,'FontSize',16);