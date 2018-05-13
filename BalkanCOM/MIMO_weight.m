%%Simulation of the paper "Power Allocation for Device-to-Device Communication
%%Underlaying Massive MIMO Multicasting Networks"
%ZF Distributed

clear all ; close all; clc


%% Parameters
R    = 200;     % Cell radius m
RD2D = 24;      % maximum D2D range
N    = 100;     % Number of base station antennas
K    = 5;       % Number of CUEs
D    = 10;       % Number of D2D pairs
Pbs  = 30;      % Base station maximum power dBm
Pbsr = 10^(Pbs/10)/1000;
Pd   = 13;      % Mobile terminal maximum power dBm
Pdr  = 10^(Pd/10)/1000;
a    = 3;       % Path loss exponen \alpha
No   = 10^(-10);
Yth  = 10.96;   % Threshold SINR for CUEs dB
Ythr = 10^(Yth/10);
weight=1;

iter=1       % Number of montecarlo iterations
iter_algorithm= 3 % iteration number for algorithm
Q=2; %iteration number for algorithm 2

%% Allocate positions of CUE and D2D pairs.
CUEratearray=[];
CUESNRarray=[];
D2Dratearray=[];
D2DSNRarray=[];
lambdaarray=[];
wn=[];
warray=[];
theta_array=[];

D2Dratearray2=[];
prcntg=[];
tic
D2Dratearray_varyingD=zeros(iter,10);
CUEratearray_varyingD=zeros(iter,10);

        %Ythr = 10^(Yth/10)*ones(K,1);

        
        
        distCUE=R*sqrt(rand(K,1)); %distance from enodeB to CUE's
        angles=2*pi*rand(K,1);
        CUEx     = distCUE.*cos(angles); % x coordinate of ith CUE
        CUEy     = distCUE.*sin(angles); % y coordinate of yth CUE
        distrxD=R*sqrt(rand(D,1)); %distance of D2D rx from enodeB
        anglerxD=2*pi*sqrt(rand(D,1));
        D_rxX    = distrxD.*cos(anglerxD); % x coordinate of ith D2D transmitter
        D_rxY    = distrxD.*sin(anglerxD); % y coordinate of ith D2D transmitter
        disttxD  = RD2D*sqrt(rand(D,1));
        %disttxD  = RD2D*rand(D,1);
        angletxD = 2*pi*rand(D,1);
        D_txX    = D_rxX + disttxD.*cos(angletxD);% x coordinate of ith D2D receiver
        D_txY    = D_rxY + disttxD.*sin(angletxD); % y coordinate of ith D2D receiver
        distCUED2D=zeros(D,K);
        for d=1:D %distance between D2D tx d and CUE k
            for k=1:K
                distCUED2D(d,k)=sqrt((CUEx(k)-D_txX(d))^2+(CUEy(k)-D_txY(d))^2);
            end
        end
        distDD=zeros(D,D);
        for d1=1:D
            for d2=1:D
                distDD(d1,d2)=sqrt((D_txX(d1)-D_rxX(d2))^2+(D_txY(d1)-D_rxY(d2))^2);
            end
        end
        %% ChannelMatrices
        h=zeros(N,K);
        betah=sqrt((distCUE).^(-a));
        for k=1:K %channel gain from enodeB to CUE's
            h(:,k)=sqrt(1/2)*(randn(N,1)+1i*randn(N,1))*betah(k);
        end
        g=zeros(D,K); %channel gain from D2D tx d to CUE k
        betag=sqrt((distCUED2D).^(-a));
        g=sqrt(1/2)*(randn(D,K)+1i*randn(D,K)).*betag;
        rho=zeros(D,D); %channel gain from D2D tx d1 to D2D rx d2
        betarho=sqrt((distDD).^(-a));
        rho=sqrt(1/2)*(randn(D,D)+1i*randn(D,D)).*betarho;
        f=zeros(N,D); %channel gain from enodeB to D2D rx d
        betaf=sqrt((distrxD).^(-a));
        for d=1:D
            f(:,d)=sqrt(1/2)*(randn(N,1)+1i*randn(N,1))*betaf(d);
        end
        %%Simple Precoders
        temporary=null(f');
        wZF=temporary(:,1);
        w=wZF;
        SNRCUENoD2D=zeros(K,1);
        for k=1:K %SNR of CUEs without any D2D transmission
            SNRCUENoD2D(k)=Pbsr*abs(h(:,k)'*w)^2/No;
        end
        
        %%Distributed
        tmp=1;
        lambdaarray=[];
        lambda=ones(D+1,1);
        for ii=-4.3:0.01:-3.8           
            weight=10^ii;
            for d=1:D
                InterfD2D_fromd(d)=sum(lambda(1:D).*Pbsr.*(abs(rho(:,d)).^2))-lambda(d)*Pdr*(abs(rho(d,d))^2);
                InterfD2D(d)=(No+lambda(D+1)*Pbsr*(abs(f(:,d)'*w))^2+InterfD2D_fromd(d));
                lambda(d)=Pdr*((abs(rho(d,d)))^2/(weight*Pd*sum(abs(g(d,:)).^2))-1)*InterfD2D(d);
                if lambda(d)>1
                    lambda(d)=1;
                elseif lambda(d)<0
                    lambda(d)=0;
                end
            end
            lambdaarray(tmp,:)=lambda
            tmp=tmp+1;
        end

        

