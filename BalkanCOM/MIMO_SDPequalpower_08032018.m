%%Simulation of the paper "Power Allocation for Device-to-Device Communication
%%Underlaying Massive MIMO Multicasting Networks"
%ZF and Equal power

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


iter=100       % Number of montecarlo iterations
iter_algorithm= 10 % iteration number for algorithm
Q=10; %iteration number for algorithm 2

% R    = 200;     % Cell radius m
% RD2D = 24;      % maximum D2D range
% N    = 10;     % Number of base station antennas
% K    = 5;       % Number of CUEs
% D    = 10;       % Number of D2D pairs
% Pbs  = 30;      % Base station maximum power dBm
% Pbsr = 10^(Pbs/10)/1000;
% Pd   = 13;      % Mobile terminal maximum power dBm
% Pdr  = 10^(Pd/10)/1000;
% a    = 3;       % Path loss exponen \alpha
% No   = 10^(-10);
% Yth  = 10.96;   % Threshold SINR for CUEs dB
% Ythr = 10^(Yth/10);
% 
% iter=100       % Number of montecarlo iterations
% iter_algorithm= 10 % iteration number for algorithm
% Q=10; %iteration number for algorithm 2

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
for monte=1:iter 
    monte
    for D=7:10
        
        Ythr = 10^(Yth/10)*ones(K,1);       
        rand('state',monte)
        randn('state',monte)
        %monte
        
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
        
        %         temporary=null(f');
        %         wZF=temporary(:,1);
        %         w=wZF;
        
        %%SDP
        lambda=ones(D+1,1); %initialize power factors of D2D tx and one enodeB
        while 1==1
            Yk=zeros(1,K);
            for k=1:K
                Yk(k)=(Ythr(1)/(Pbsr*lambda(D+1)))*(sum(lambda(1:D).*Pdr.*(abs(g(:,k)).^2))+No);
                %Yk(k)=0;
            end
            
            check=10;
            while check>D
                F{check}=0;
                check=check-1;
            end
            for d=1:D
                F{d}=f(:,d)*f(:,d)';
            end
        
        
            lambda'
            cvx_begin sdp quiet   
            variable X(N,N) hermitian
            cost=0;
            for d=1:D
                A(d)=(1/(Pbsr*lambda(D+1)))*(sum(lambda(1:D).*Pdr.*(abs(rho(:,d)).^2))+No);
                B(d)=(1/(Pbsr*lambda(D+1)))*(sum(lambda(1:D).*Pdr.*(abs(rho(:,d)).^2))-lambda(d)*Pdr*(abs(rho(d,d))^2)+No);                
                cost= cost+((A(d)-B(d))/(A(d)*B(d)))*real(trace((f(:,d)*f(:,d)')*X));
            end
            minimize cost
            trace(X) == 1;
            X >= 0;
            
            for k=1:K
                real(trace((h(:,k)*h(:,k)')*X))-Yk(k)>= 0
            end
            cvx_end
            
            fprintf(1,'Done! \n');
            opt2 = cvx_optval;
            if opt2 ~= Inf
                break
            else
                lambda(1:D)=max(0,lambda(1:D)-0.05);
            end
        
        end
        if opt2 ~= Inf
            
            [Vx,Dx]=eig(X);
            Vxx=eig(X)
            E=exp(i*rand(rank(Vx),N)*2*pi);
            W=Vx*Dx^(1/2)*E
            for z=1:size(Vx,2)
                for d=1:D %calculate the d2D rates before precoding
                    wx=W(:,z);
                    InterfD2D(d)=sum(lambda(1:D).*Pdr.*(abs(rho(:,d)).^2))-lambda(d)*Pdr*(abs(rho(d,d))^2);
                    SNRD2D(d)=lambda(d)*Pdr*(abs(rho(d,d)))^2/(No+lambda(D+1)*Pbsr*(abs(f(:,d)'*wx))^2+InterfD2D(d));
                    RateD2D(d)=log2(1+SNRD2D(d));
                end
                sumRateD2Dx(z)=sum(RateD2D);
            end
            optimumprecoderindex=find(sumRateD2Dx==max(sumRateD2Dx));
            w=W(:,optimumprecoderindex(1))
            
            
            
            %%%%%%%
            
            SNRCUENoD2D=zeros(K,1);
            for k=1:K %SNR of CUEs without any D2D transmission
                SNRCUENoD2D(k)=Pbsr*abs(h(:,k)'*w)^2/No;
            end
            

            
            for k=1:K %SNR of CUEs without any D2D transmission
                InterfCUE(k)=sum(Pdr*lambda(1:D).*(abs(g(:,k)).^2));
                SNRCUE(k)= lambda(D+1)*Pbsr*(abs(h(:,k)'*w)^2)/(InterfCUE(k)+No);
                RateCUE(k)=log2(1+SNRCUE(k));
            end

            SNRD2D=zeros(D,1);
            InterfD2D=zeros(D,1);
            RateD2D=zeros(D,1);
            
            for d=1:D %calculate the d2D rates before precoding
                InterfD2D(d)=sum(lambda(1:D).*Pdr.*(abs(rho(:,d)).^2))-lambda(d)*Pdr*(abs(rho(d,d))^2);
                SNRD2D(d)=lambda(d)*Pdr*(abs(rho(d,d)))^2/(No+lambda(D+1)*Pbsr*(abs(f(:,d)'*w))^2+InterfD2D(d));
                RateD2D(d)=log2(1+SNRD2D(d));
            end
            %RateD2D
            D2Dratearray_varyingD(monte,D)=sum(RateD2D);
            CUEratearray_varyingD(monte,D)=sum(RateCUE);
                        
            CUESNRarray=[CUESNRarray;SNRCUE];
        else
            disp(['Inf! for']);
            monte
        end
    end
        
end

for d=1:D
    totalNonzeroresults(d)=sum(D2Dratearray_varyingD(:,d)~=0);
end
d2dresult=sum(D2Dratearray_varyingD)./totalNonzeroresults
meanpercent=mean(prcntg)
ZF_EqualPower_D2Dsumrate_withoutprecoder=d2dresult %mean(D2Dratearray_varyingD); %average D2Dsum rate with BF precoder
ZF_CUEsumrate_withoutprecoder=mean(CUEratearray_varyingD); %average sum rate with BF precoder

%Plot figures
figure('Name','Sum rate of D2D pairs','NumberTitle','off')
plot(ZF_EqualPower_D2Dsumrate_withoutprecoder,'-o','color','g')
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
legend('Prec=ZF PA=Uniform')

%figure
%cdfplot(10*log10(CUESNRarray))

toc
