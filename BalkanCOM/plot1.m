qiteration=[0.5141    0.1210    0.0243    0.0113    0.0056    0.0035    0.0022    0.0018    0.0012;
0.3678    0.1194    0.0311    0.0170    0.0087    0.0071    0.0044    0.0039    0.0024;
0.3191    0.1016    0.0336    0.0218    0.0152    0.0110    0.0081    0.0058    0.0052]

figure('Name','GD algorithm (Q)','NumberTitle','off')
hold on
plot(qiteration(1,:),'-o','color','g')
plot(qiteration(2,:),'-x','color','r')
plot(qiteration(3,:),'-O','color','b')
grid;
xlabel('Iterations');
ylabel('Increase in D2D Rates');
title('GD Algorithm iterations Q')
legend('D=3','D=5','D=10')

sumrate=[3.6069    6.3761    8.3038    9.7449   10.7351   11.6572   12.3428   12.8623   12.5605   12.5340;
8.3038   12.5725   14.0934   14.4353   14.5988   14.6801   14.7316   14.7637   14.7902   14.8074;
11.3328   14.0698   16.1280   17.8363   18.3573   19.1946   18.9876   20.4497   20.2047   19.8549]


%SDP precoder [equal,heuristic,distirbuted]
sdp_sumrate=[14.98547926	20.00509608	24.65018728	27.35989908	29.84853938	32.76093084	35.4009039	39.39821821	39.68643156	41.40361736;
[14.3042371254638 18.4633161055943 23.5873679637112 27.8286974717192 30.3094679937993 33.1915606711982 35.4806659358300 40.1670133734077 42.7063694892700 44.6040059927449];
[14.5885682315325 19.5475189867560 24.3314833603472 26.3468323106870 28.8806466520951 32.5881597414833 34.5751575868490 37.3911627842767 39.0675913077247 41.3610599666360]];
%    
figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sdp_sumrate(1,:),'-o','color','g')
plot(sdp_sumrate(2,:),'-x','color','r')
plot(sdp_sumrate(3,:),'-o','color','b')

legend('Prec=SDP PA=EqualPower','Prec=SDP PA=Heuristic','Prec=SDP PA=Distributed')



%equalpower [BF,GD,ZF,SDP]
sumrate=[3.6069    6.3761    8.3038    9.7449   10.7351   11.6572   12.3428   12.8623   12.5605   12.5340;
9.3105   10.7724   12.9295   13.4518   14.5265   16.0918   17.7414   18.7408   18.3661   17.7502;
15.0566   19.3721   12.5589   11.9142   11.5111   10.8016    9.3173    8.7468    7.3884    5.9136;
14.98547926	20.00509608	24.65018728	27.35989908	29.84853938	32.76093084	35.4009039	39.39821821	39.68643156	41.40361736
];

figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumrate(1,:),'-o','color','g')
plot(sumrate(2,:),'-x','color','r')
plot(sumrate(3,:),'-o','color','b')
plot(sumrate(4,:),'-o','color','black')
legend('Prec=BF PA=EqualPower','Prec=GD PA=EqualPower','Prec=ZF PA=EqualPower', 'Prec=SDP PA=EqualPower')



%heuristic PA [BF,GD,ZF,SDP]
sumrate=[3.6884    6.3970    8.4408   10.0573   11.1432   12.4611   13.3766   14.3554   14.2104   14.8951;
11.0141   13.6273   16.1718   18.0902   19.8940   21.1184   24.4496   24.5472   26.1927   27.6691;
12.8247   16.3828   20.0727   23.5695   26.5800   29.4912   31.6229   34.6557   35.2572   36.9664;
[14.3042371254638 18.4633161055943 23.5873679637112 27.8286974717192 30.3094679937993 33.1915606711982 35.4806659358300 40.1670133734077 42.7063694892700 44.6040059927449]];
%[15.8704454501286 18.9081573268708 25.5113397612820 27.4279564196706 28.6670219416573 31.9793208426105 34.8017708343136 37.0314168675697 39.4096141891599 39.1650308420705]];


figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumrate(1,:),'-o','color','g')
plot(sumrate(2,:),'-x','color','r')
plot(sumrate(3,:),'-o','color','b')
plot(sumrate(4,:),'-o','color','black')
legend('Prec=BF PA=Heuristic','Prec=GD PA=Heuristic','Prec=ZF PA=Heuristic', 'Prec=SDP PA=Heuristic')



%distributed [BF,GD,ZF]
% sumrate=[2.8958    5.6318    8.4387   10.8673   13.2251   15.5939   17.9948   20.0941   21.4168   23.8128;
% 8.3562   12.6537   15.6330   17.1841   19.8122   23.0322   27.1281   31.3520   31.5267   35.2771;
% 11.5370   17.8900   22.4434   26.0167   29.3324   32.1792   35.2823   37.5380   39.1853   41.7296];
sumrate=[3.8434    7.4503    9.4073   11.1858   13.7501   16.1436   17.6338   19.1978   21.1137   22.8800;
 11.3357   13.6669   17.2661   20.4201   23.0011   25.8347   29.9705   31.0433   31.3365   31.9173;
12.8245   17.0292   20.1790   22.9249   25.0765   27.1900   28.5284   30.4852   31.4478   32.6499;
[14.5885682315325 19.5475189867560 24.3314833603472 26.3468323106870 28.8806466520951 32.5881597414833 34.5751575868490 37.3911627842767 39.0675913077247 41.3610599666360]];
%[14.7944758182132 15.2184286764504 24.0202708502226 27.6594021343844 29.0215987894158 29.1509367120521 31.7845090741225 35.1118943220703 41.9046127898339 46.7697243840506]];


figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumrate(1,:),'-o','color','g')
plot(sumrate(2,:),'-x','color','r')
plot(sumrate(3,:),'-o','color','blue')
plot(sumrate(4,:),'-o','color','black')

legend('Prec=BF PA=Distributed','Prec=GD PA=Distributed','Prec=ZF PA=Distributed', 'Prec=SDP PA=Distributed')

%GD 
sumrate=[8.3038   12.5725   14.0934   14.4353   14.5988   14.6801   14.7316   14.7637   14.7902   14.8074;
11.0141   13.6273   16.1718   18.0902   19.8940   21.1184   24.4496   24.5472   26.1927   27.6691;
11.3357   13.6669   17.2661   20.4201   23.0011   25.8347   29.9705   31.0433   31.3365   31.9173];

figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumrate(1,:),'-o','color','g')
plot(sumrate(2,:),'-x','color','r')
plot(sumrate(3,:),'-o','color','b')

legend('Prec=GD PA=Uniform','Prec=GD PA=Heuristic','Prec=GD PA=Distributed')




figure
hold on
cdfplot(10*log10(CUESNRarrayZF))
cdfplot(10*log10(CUESNRarrayBF))
cdfplot(10*log10(CUESNRarrayGD))
legend('Prec=ZF PA=Uniform','Prec=BF PA=Uniform','Prec=GD PA=Uniform','Prec=SDP PA=Uniform')


%%

%sumrate new equal power [ BF,GD,ZF,SDP] N=10

sumratenew=[2.9903    4.4446    5.3944    5.6360    5.6978    5.7684    5.1446    5.0058    4.5361    4.5735;
5.4484    6.9438    7.9960    8.3860    8.5528    8.7246    7.6919    7.5275    6.6471    6.4847;
11.0729   12.2501   12.7590   12.3692   11.8001   11.3074    9.0209    8.9716    7.4522   5.9136;
2.9241    3.6204    3.7586    2.8834    2.5983    2.5020    1.8212    1.6019    2.3572    1.9933]
figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumratenew(1,:),'-o','color','g')
plot(sumratenew(2,:),'-x','color','r')
plot(sumratenew(3,:),'-o','color','b')
plot(sumratenew(4,:),'-o','color','black')
legend('Prec=BF PA=Uniform','Prec=GD PA=Uniform','Prec=ZF PA=Uniform', 'Prec=SDP PA=Uniform')

sumratenew=[3.4253    5.7797    7.8866    9.5364   11.1284   12.4319   13.4676   14.0613   15.3246   16.7529;
5.9938    8.3761   10.6898   12.5165   14.3504   15.6355   16.8313   17.4247   18.2613   19.8309;
13.8758   18.6740   23.2829   26.4300   31.4669   34.2649   36.8200   38.8381   41.8664  43.1650
13.7929   17.8346   23.6389   26.6982   31.3047   33.3598   39.5917   39.4059   41.3551   44.2960]

figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Nonzero Lambda average Sum rate of D2D pairs (bits/sec) ');
title('Sum rate of D2D pairs')
plot(sumratenew(1,:),'-o','color','g')
plot(sumratenew(2,:),'-x','color','r')
plot(sumratenew(3,:),'-o','color','b')
plot(sumratenew(4,:),'-o','color','black')
legend('Prec=BF PA=Uniform','Prec=GD PA=Uniform','Prec=ZF PA=Uniform', 'Prec=SDP PA=Uniform')


%BF 10.04.2018
BFvsPA=[3.6069    6.3761    8.3038    9.7449   10.7351   11.6572   12.3428   12.8623   12.5605   12.5340;
    3.6884    6.3970    8.4408   10.0573   11.1432   12.4611   13.3766   14.3554   14.2104   14.8951;
    [3.68839493927662 6.83899613896662 9.53960260141339 11.5928490236660 13.4238851175011 15.3005573454332 17.1300186926143 18.7693943273552 19.2842380102654 20.9668080461721]];
figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(BFvsPA(1,:),'-o','color','g')
plot(BFvsPA(2,:),'-x','color','r')
plot(BFvsPA(3,:),'-o','color','b')

legend('Prec=BF PA=Uniform','Prec=BF PA=GD','Prec=BF PA=Distributed')

%GD
BFvsPA=[9.3105   10.7724   12.9295   13.4518   14.5265   16.0918   17.7414   18.7408   18.3661   17.7502;
    11.0141   13.6273   16.1718   18.0902   19.8940   21.1184   24.4496   24.5472   26.1927   27.6691;
    11.3357   13.6669   17.2661   20.4201   23.0011   25.8347   29.9705   31.0433   31.3365   31.9173;]

figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumrate(1,:),'-o','color','g')
plot(sumrate(2,:),'-x','color','r')
plot(sumrate(3,:),'-o','color','b')

legend('Prec=GD PA=Uniform','Prec=GD PA=Heuristic','Prec=GD PA=Distributed')


%ZF
sumrate=[15.0566   19.3721   12.5589   11.9142   11.5111   10.8016    9.3173    8.7468    7.3884    5.9136;
    12.8247   16.3828   20.0727   23.5695   26.5800   29.4912   31.6229   34.6557   35.2572   36.9664;
    12.8245   17.0292   20.1790   22.9249   25.0765   27.1900   28.5284   30.4852   31.4478   32.6499;];

figure('Name','Sum rate of D2D pairs','NumberTitle','off')
hold on
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
plot(sumrate(1,:),'-o','color','g')
plot(sumrate(2,:),'-x','color','r')
plot(sumrate(3,:),'-o','color','b')

legend('Prec=ZF PA=Uniform','Prec=ZF PA=Heuristic','Prec=ZF PA=Distributed')
