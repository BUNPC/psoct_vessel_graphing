load matlab;

%%
s = EllipseFit3DConstrained_dab( V,72,77,40,1);

%%

ShowLocalDataWithSE(V,s,'foo')