%% colormaps
function [Vcmap, Vzcmap, Dcmap, Mfcmap, Rcmap, g1OCTAcmap]=Colormaps_DLSOCT
% h=gca; custColormap=colormap(h);
Vcmap=[0         0         0
         0         0    0.4375
         0         0    0.8750
         0         0    0.9000
         0         0    0.9250
         0         0    0.9500
         0         0    0.9750
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0
    0.9375         0         0
    0.8750         0         0
    0.8125         0         0
    0.7500         0         0
    0.6875         0         0
    0.6250         0         0
    0.5625         0         0
    0.5000         0         0];
% Vcmap=[0         0         0
%     0.0176         0         0
%     0.0351         0         0
%     0.0527         0         0
%     0.0703         0         0
%     0.0878         0         0
%     0.1054         0         0
%     0.1230         0         0
%     0.1405         0         0
%     0.1581         0         0
%     0.1757         0         0
%     0.1932         0         0
%     0.2108         0         0
%     0.2283         0         0
%     0.2351         0         0
%     0.2419         0         0
%     0.2487         0         0
%     0.2554         0         0
%     0.2622         0         0
%     0.2690         0         0
%     0.2757         0         0
%     0.2825         0         0
%     0.2893         0         0
%     0.2960         0         0
%     0.3028         0         0
%     0.3096         0         0
%     0.3163         0         0
%     0.3231         0         0
%     0.3299         0         0
%     0.3366         0         0
%     0.3434         0         0
%     0.3502         0         0
%     0.3570         0         0
%     0.3637         0         0
%     0.3705         0         0
%     0.3773         0         0
%     0.3840         0         0
%     0.3908         0         0
%     0.3976         0         0
%     0.4043         0         0
%     0.4111         0         0
%     0.4179         0         0
%     0.4246         0         0
%     0.4314         0         0
%     0.4382         0         0
%     0.4450         0         0
%     0.4517         0         0
%     0.4585         0         0
%     0.4653         0         0
%     0.4720         0         0
%     0.4788         0         0
%     0.4856         0         0
%     0.4923         0         0
%     0.4991         0         0
%     0.5059         0         0
%     0.5126         0         0
%     0.5194         0         0
%     0.5262         0         0
%     0.5329         0         0
%     0.5397         0         0
%     0.5465         0         0
%     0.5533         0         0
%     0.5600         0         0
%     0.5668         0         0
%     0.5736         0         0
%     0.5803         0         0
%     0.5871         0         0
%     0.5939         0         0
%     0.6006         0         0
%     0.6074         0         0
%     0.6142         0         0
%     0.6209         0         0
%     0.6277         0         0
%     0.6345         0         0
%     0.6412         0         0
%     0.6480         0         0
%     0.6548         0         0
%     0.6616         0         0
%     0.6683         0         0
%     0.6751         0         0
%     0.6819         0         0
%     0.6886         0         0
%     0.6954         0         0
%     0.7022         0         0
%     0.7089         0         0
%     0.7157         0         0
%     0.7225         0         0
%     0.7292         0         0
%     0.7360         0         0
%     0.7428         0         0
%     0.7496         0         0
%     0.7563         0         0
%     0.7631         0         0
%     0.7699         0         0
%     0.7766         0         0
%     0.7834         0         0
%     0.7902         0         0
%     0.7969         0         0
%     0.8037         0         0
%     0.8105         0         0
%     0.8172         0         0
%     0.8240         0         0
%     0.8308         0         0
%     0.8375         0         0
%     0.8443         0         0
%     0.8511         0         0
%     0.8579         0         0
%     0.8646         0         0
%     0.8714         0         0
%     0.8782         0         0
%     0.8849         0         0
%     0.8917         0         0
%     0.8985         0         0
%     0.9052         0         0
%     0.9120         0         0
%     0.9188         0         0
%     0.9255         0         0
%     0.9323         0         0
%     0.9391         0         0
%     0.9458         0         0
%     0.9526         0         0
%     0.9594         0         0
%     0.9662         0         0
%     0.9729         0         0
%     0.9797         0         0
%     0.9865         0         0
%     0.9932         0         0
%     1.0000         0         0];
Vzcmap=[0.6784    0.9216    1.0000
    0.6262    0.9276    0.9231
    0.5741    0.9336    0.8462
    0.5219    0.9397    0.7692
    0.4697    0.9457    0.6923
    0.4175    0.9517    0.6154
    0.3653    0.9578    0.5385
    0.3131    0.9638    0.4615
    0.2609    0.9698    0.3846
    0.2087    0.9759    0.3077
    0.1566    0.9819    0.2308
    0.1044    0.9879    0.1538
    0.0522    0.9940    0.0769
         0    1.0000         0
         0    0.9442         0
         0    0.8885         0
         0    0.8327         0
         0    0.7769         0
         0    0.7211         0
         0    0.6654         0
         0    0.6096         0
         0    0.5538         0
         0    0.4980         0
         0    0.4482         0
         0    0.3984         0
         0    0.3486         0
         0    0.2988         0
         0    0.2490         0
         0    0.1992         0
         0    0.1494         0
         0    0.0996         0
         0    0.0498         0
         0         0         0
    0.0909         0         0
    0.1818         0         0
    0.2727         0         0
    0.3636         0         0
    0.4545         0         0
    0.5455         0         0
    0.6364         0         0
    0.7273         0         0
    0.8182         0         0
    0.9091         0         0
    1.0000         0         0
    1.0000    0.1000         0
    1.0000    0.2000         0
    1.0000    0.3000         0
    1.0000    0.4000         0
    1.0000    0.5000         0
    1.0000    0.6000         0
    1.0000    0.7000         0
    1.0000    0.8000         0
    1.0000    0.9000         0
    1.0000    1.0000         0
    0.9973    0.9973    0.0973
    0.9945    0.9945    0.1945
    0.9918    0.9918    0.2918
    0.9890    0.9890    0.3890
    0.9863    0.9863    0.4863
    0.9835    0.9835    0.5835
    0.9808    0.9808    0.6808
    0.9780    0.9780    0.7780
    0.9753    0.9753    0.8753
    0.9725    0.9725    0.9725];
% Vzcmap=[0    1.0000    1.0000
%          0    0.9893    0.9893
%          0    0.9785    0.9785
%          0    0.9678    0.9678
%          0    0.9570    0.9570
%          0    0.9463    0.9463
%          0    0.9355    0.9355
%          0    0.9248    0.9248
%          0    0.9140    0.9140
%          0    0.9033    0.9033
%          0    0.8925    0.8925
%          0    0.8818    0.8818
%          0    0.8710    0.8710
%          0    0.8603    0.8603
%          0    0.8495    0.8495
%          0    0.8388    0.8388
%          0    0.8280    0.8280
%          0    0.8173    0.8173
%          0    0.8065    0.8065
%          0    0.7958    0.7958
%          0    0.7850    0.7850
%          0    0.7743    0.7743
%          0    0.7635    0.7635
%          0    0.7528    0.7528
%          0    0.7421    0.7421
%          0    0.7313    0.7313
%          0    0.7206    0.7206
%          0    0.7098    0.7098
%          0    0.6991    0.6991
%          0    0.6883    0.6883
%          0    0.6776    0.6776
%          0    0.6668    0.6668
%          0    0.6561    0.6561
%          0    0.6453    0.6453
%          0    0.6346    0.6346
%          0    0.6238    0.6238
%          0    0.6131    0.6131
%          0    0.6023    0.6023
%          0    0.5916    0.5916
%          0    0.5808    0.5808
%          0    0.5701    0.5701
%          0    0.5593    0.5593
%          0    0.5486    0.5486
%          0    0.5378    0.5378
%          0    0.5271    0.5271
%          0    0.5163    0.5163
%          0    0.5056    0.5056
%          0    0.4948    0.4948
%          0    0.4841    0.4841
%          0    0.4734    0.4734
%          0    0.4626    0.4626
%          0    0.4519    0.4519
%          0    0.4411    0.4411
%          0    0.4304    0.4304
%          0    0.4196    0.4196
%          0    0.4089    0.4089
%          0    0.3981    0.3981
%          0    0.3874    0.3874
%          0    0.3766    0.3766
%          0    0.3659    0.3659
%          0    0.3551    0.3551
%          0    0.3444    0.3444
%          0    0.3336    0.3336
%          0    0.3229    0.3229
%          0    0.3121    0.3121
%          0    0.3014    0.3014
%          0    0.2906    0.2906
%          0    0.2799    0.2799
%          0    0.2691    0.2691
%          0    0.2584    0.2584
%          0    0.2476    0.2476
%          0    0.2369    0.2369
%          0    0.2262    0.2262
%          0    0.2154    0.2154
%          0    0.2047    0.2047
%          0    0.1939    0.1939
%          0    0.1832    0.1832
%          0    0.1724    0.1724
%          0    0.1552    0.1552
%          0    0.1379    0.1379
%          0    0.1207    0.1207
%          0    0.1034    0.1034
%          0    0.0862    0.0862
%          0    0.0690    0.0690
%          0    0.0517    0.0517
%          0    0.0345    0.0345
%          0    0.0172    0.0172
%          0         0         0
%     0.0281         0         0
%     0.0562         0         0
%     0.0844         0         0
%     0.1125         0         0
%     0.1406         0         0
%     0.1687         0         0
%     0.1969         0         0
%     0.2250         0         0
%     0.2492         0         0
%     0.2734         0         0
%     0.2977         0         0
%     0.3219         0         0
%     0.3461         0         0
%     0.3703         0         0
%     0.3945         0         0
%     0.4187         0         0
%     0.4430         0         0
%     0.4672         0         0
%     0.4914         0         0
%     0.5156         0         0
%     0.5398         0         0
%     0.5641         0         0
%     0.5883         0         0
%     0.6125         0         0
%     0.6367         0         0
%     0.6609         0         0
%     0.6852         0         0
%     0.7094         0         0
%     0.7336         0         0
%     0.7578         0         0
%     0.7820         0         0
%     0.8062         0         0
%     0.8305         0         0
%     0.8547         0         0
%     0.8789         0         0
%     0.9031         0         0
%     0.9273         0         0
%     0.9516         0         0
%     0.9758         0         0
%     1.0000         0         0];
Dcmap=[       0         0         0
    0.0024    0.0024         0
    0.0048    0.0048         0
    0.0072    0.0072         0
    0.0096    0.0096         0
    0.0120    0.0120         0
    0.0144    0.0144         0
    0.0168    0.0168         0
    0.0192    0.0192         0
    0.0216    0.0216         0
    0.0240    0.0240         0
    0.0264    0.0264         0
    0.0288    0.0288         0
    0.0312    0.0312         0
    0.0336    0.0336         0
    0.0359    0.0359         0
    0.0383    0.0383         0
    0.0407    0.0407         0
    0.0431    0.0431         0
    0.0455    0.0455         0
    0.0479    0.0479         0
    0.0503    0.0503         0
    0.0527    0.0527         0
    0.0551    0.0551         0
    0.0642    0.0642         0
    0.0733    0.0733         0
    0.0824    0.0824         0
    0.0915    0.0915         0
    0.1005    0.1005         0
    0.1096    0.1096         0
    0.1187    0.1187         0
    0.1278    0.1278         0
    0.1369    0.1369         0
    0.1460    0.1460         0
    0.1551    0.1551         0
    0.1641    0.1641         0
    0.1732    0.1732         0
    0.1823    0.1823         0
    0.1914    0.1914         0
    0.2005    0.2005         0
    0.2096    0.2096         0
    0.2187    0.2187         0
    0.2277    0.2277         0
    0.2368    0.2368         0
    0.2459    0.2459         0
    0.2550    0.2550         0
    0.2641    0.2641         0
    0.2732    0.2732         0
    0.2823    0.2823         0
    0.2913    0.2913         0
    0.3004    0.3004         0
    0.3095    0.3095         0
    0.3186    0.3186         0
    0.3277    0.3277         0
    0.3368    0.3368         0
    0.3459    0.3459         0
    0.3549    0.3549         0
    0.3640    0.3640         0
    0.3731    0.3731         0
    0.3822    0.3822         0
    0.3913    0.3913         0
    0.4004    0.4004         0
    0.4094    0.4094         0
    0.4185    0.4185         0
    0.4276    0.4276         0
    0.4367    0.4367         0
    0.4458    0.4458         0
    0.4549    0.4549         0
    0.4640    0.4640         0
    0.4730    0.4730         0
    0.4821    0.4821         0
    0.4912    0.4912         0
    0.5003    0.5003         0
    0.5094    0.5094         0
    0.5185    0.5185         0
    0.5276    0.5276         0
    0.5366    0.5366         0
    0.5457    0.5457         0
    0.5548    0.5548         0
    0.5639    0.5639         0
    0.5730    0.5730         0
    0.5821    0.5821         0
    0.5912    0.5912         0
    0.6002    0.6002         0
    0.6093    0.6093         0
    0.6184    0.6184         0
    0.6275    0.6275         0
    0.6366    0.6366         0
    0.6457    0.6457         0
    0.6548    0.6548         0
    0.6638    0.6638         0
    0.6729    0.6729         0
    0.6820    0.6820         0
    0.6911    0.6911         0
    0.7002    0.7002         0
    0.7093    0.7093         0
    0.7184    0.7184         0
    0.7274    0.7274         0
    0.7365    0.7365         0
    0.7456    0.7456         0
    0.7547    0.7547         0
    0.7638    0.7638         0
    0.7729    0.7729         0
    0.7820    0.7820         0
    0.7910    0.7910         0
    0.8001    0.8001         0
    0.8092    0.8092         0
    0.8183    0.8183         0
    0.8274    0.8274         0
    0.8365    0.8365         0
    0.8455    0.8455         0
    0.8546    0.8546         0
    0.8637    0.8637         0
    0.8728    0.8728         0
    0.8819    0.8819         0
    0.8910    0.8910         0
    0.9001    0.9001         0
    0.9091    0.9091         0
    0.9182    0.9182         0
    0.9273    0.9273         0
    0.9364    0.9364         0
    0.9455    0.9455         0
    0.9546    0.9546         0
    0.9637    0.9637         0
    0.9727    0.9727         0
    0.9818    0.9818         0
    0.9909    0.9909         0
    1.0000    1.0000         0];
Mfcmap=[0         0         0
    0.0127         0         0
    0.0254         0         0
    0.0381         0         0
    0.0889         0         0
    0.1397         0         0
    0.1905         0         0
    0.2189         0         0
    0.2473         0         0
    0.2757         0         0
    0.3041         0         0
    0.3172         0         0
    0.3304         0         0
    0.3435         0         0
    0.3566         0         0
    0.3697         0         0
    0.3829         0         0
    0.3960         0         0
    0.4091         0         0
    0.4223         0         0
    0.4354         0         0
    0.4485         0         0
    0.4617         0         0
    0.4748         0         0
    0.4879         0         0
    0.5010         0         0
    0.5142         0         0
    0.5273         0         0
    0.5404         0         0
    0.5536         0         0
    0.5667         0         0
    0.5798         0         0
    0.5930         0         0
    0.6061         0         0
    0.6192         0         0
    0.6324         0         0
    0.6455         0         0
    0.6586         0         0
    0.6717         0         0
    0.6849         0         0
    0.6980         0         0
    0.7111         0         0
    0.7243         0         0
    0.7374         0         0
    0.7505         0         0
    0.7637         0         0
    0.7768         0         0
    0.7899         0         0
    0.8030         0         0
    0.8162         0         0
    0.8293         0         0
    0.8424         0         0
    0.8556         0         0
    0.8687         0         0
    0.8818         0         0
    0.8950         0         0
    0.9081         0         0
    0.9212         0         0
    0.9343         0         0
    0.9475         0         0
    0.9606         0         0
    0.9737         0         0
    0.9869         0         0
    1.0000         0         0];
Rcmap=[0    1.0000    1.0000
    0.0140    0.9860    1.0000
    0.0280    0.9720    1.0000
    0.0420    0.9580    1.0000
    0.0560    0.9440    1.0000
    0.0700    0.9300    1.0000
    0.0840    0.9160    1.0000
    0.0980    0.9020    1.0000
    0.1120    0.8880    1.0000
    0.1260    0.8740    1.0000
    0.1400    0.8600    1.0000
    0.1540    0.8460    1.0000
    0.1680    0.8320    1.0000
    0.1820    0.8180    1.0000
    0.1960    0.8040    1.0000
    0.2100    0.7900    1.0000
    0.2240    0.7760    1.0000
    0.2380    0.7620    1.0000
    0.2520    0.7480    1.0000
    0.2588    0.7412    1.0000
    0.2657    0.7343    1.0000
    0.2726    0.7274    1.0000
    0.2794    0.7206    1.0000
    0.2863    0.7137    1.0000
    0.2931    0.7069    1.0000
    0.3000    0.7000    1.0000
    0.3069    0.6931    1.0000
    0.3137    0.6863    1.0000
    0.3206    0.6794    1.0000
    0.3275    0.6725    1.0000
    0.3343    0.6657    1.0000
    0.3412    0.6588    1.0000
    0.3480    0.6520    1.0000
    0.3549    0.6451    1.0000
    0.3618    0.6382    1.0000
    0.3686    0.6314    1.0000
    0.3755    0.6245    1.0000
    0.3824    0.6176    1.0000
    0.3892    0.6108    1.0000
    0.3961    0.6039    1.0000
    0.4029    0.5971    1.0000
    0.4098    0.5902    1.0000
    0.4167    0.5833    1.0000
    0.4235    0.5765    1.0000
    0.4304    0.5696    1.0000
    0.4373    0.5627    1.0000
    0.4441    0.5559    1.0000
    0.4510    0.5490    1.0000
    0.4578    0.5422    1.0000
    0.4647    0.5353    1.0000
    0.4716    0.5284    1.0000
    0.4784    0.5216    1.0000
    0.4853    0.5147    1.0000
    0.4922    0.5078    1.0000
    0.4990    0.5010    1.0000
    0.5059    0.4941    1.0000
    0.5128    0.4872    1.0000
    0.5196    0.4804    1.0000
    0.5265    0.4735    1.0000
    0.5333    0.4667    1.0000
    0.5402    0.4598    1.0000
    0.5471    0.4529    1.0000
    0.5539    0.4461    1.0000
    0.5608    0.4392    1.0000
    0.5677    0.4323    1.0000
    0.5745    0.4255    1.0000
    0.5814    0.4186    1.0000
    0.5882    0.4118    1.0000
    0.5951    0.4049    1.0000
    0.6020    0.3980    1.0000
    0.6088    0.3912    1.0000
    0.6157    0.3843    1.0000
    0.6226    0.3774    1.0000
    0.6294    0.3706    1.0000
    0.6363    0.3637    1.0000
    0.6431    0.3569    1.0000
    0.6500    0.3500    1.0000
    0.6569    0.3431    1.0000
    0.6637    0.3363    1.0000
    0.6706    0.3294    1.0000
    0.6775    0.3225    1.0000
    0.6843    0.3157    1.0000
    0.6912    0.3088    1.0000
    0.6980    0.3020    1.0000
    0.7049    0.2951    1.0000
    0.7118    0.2882    1.0000
    0.7186    0.2814    1.0000
    0.7255    0.2745    1.0000
    0.7324    0.2676    1.0000
    0.7392    0.2608    1.0000
    0.7461    0.2539    1.0000
    0.7529    0.2471    1.0000
    0.7598    0.2402    1.0000
    0.7667    0.2333    1.0000
    0.7735    0.2265    1.0000
    0.7804    0.2196    1.0000
    0.7873    0.2127    1.0000
    0.7941    0.2059    1.0000
    0.8010    0.1990    1.0000
    0.8078    0.1922    1.0000
    0.8147    0.1853    1.0000
    0.8216    0.1784    1.0000
    0.8284    0.1716    1.0000
    0.8353    0.1647    1.0000
    0.8422    0.1578    1.0000
    0.8490    0.1510    1.0000
    0.8559    0.1441    1.0000
    0.8627    0.1373    1.0000
    0.8696    0.1304    1.0000
    0.8765    0.1235    1.0000
    0.8833    0.1167    1.0000
    0.8902    0.1098    1.0000
    0.8971    0.1029    1.0000
    0.9039    0.0961    1.0000
    0.9108    0.0892    1.0000
    0.9176    0.0824    1.0000
    0.9245    0.0755    1.0000
    0.9314    0.0686    1.0000
    0.9382    0.0618    1.0000
    0.9451    0.0549    1.0000
    0.9520    0.0480    1.0000
    0.9588    0.0412    1.0000
    0.9657    0.0343    1.0000
    0.9725    0.0275    1.0000
    0.9794    0.0206    1.0000
    0.9863    0.0137    1.0000
    0.9931    0.0069    1.0000
    1.0000         0    1.0000];

g1OCTAcmap=[       0    1.0000    1.0000
    0.0016    0.9750    0.9750
    0.0032    0.9500    0.9500
    0.0048    0.9250    0.9250
    0.0064    0.9000    0.9000
    0.0080    0.8750    0.8750
    0.0096    0.8500    0.8500
    0.0112    0.8250    0.8250
    0.0128    0.8000    0.8000
    0.0144    0.7750    0.7750
    0.0160    0.7500    0.7500
    0.0176    0.7250    0.7250
    0.0192    0.7000    0.7000
    0.0208    0.6750    0.6750
    0.0224    0.6500    0.6500
    0.0240    0.6250    0.6250
    0.0257    0.5983    0.5983
    0.0274    0.5715    0.5715
    0.0291    0.5448    0.5448
    0.0308    0.5180    0.5180
    0.0326    0.4913    0.4913
    0.0343    0.4646    0.4646
    0.0360    0.4378    0.4378
    0.0377    0.4111    0.4111
    0.0394    0.3844    0.3844
    0.0411    0.3576    0.3576
    0.0428    0.3309    0.3309
    0.0464    0.2757    0.2757
    0.0499    0.2206    0.2206
    0.0534    0.1654    0.1654
    0.0569    0.1103    0.1103
    0.0605    0.0551    0.0551
    0.0640         0         0
    0.0820         0         0
    0.1000         0         0
    0.1180         0         0
    0.1360         0         0
    0.2440         0         0
    0.3520         0         0
    0.4600         0         0
    0.5680         0         0
    0.6760         0         0
    0.7840         0         0
    0.8920         0         0
    1.0000         0         0
    1.0000    0.0833         0
    1.0000    0.1667         0
    1.0000    0.2500         0
    1.0000    0.3333         0
    1.0000    0.4167         0
    1.0000    0.5000         0
    1.0000    0.5833         0
    1.0000    0.6667         0
    1.0000    0.7500         0
    1.0000    0.8333         0
    1.0000    0.9167         0
    1.0000    1.0000         0
    1.0000    1.0000    0.1429
    1.0000    1.0000    0.2857
    1.0000    1.0000    0.4286
    1.0000    1.0000    0.5714
    1.0000    1.0000    0.7143
    1.0000    1.0000    0.8571
    1.0000    1.0000    1.0000];