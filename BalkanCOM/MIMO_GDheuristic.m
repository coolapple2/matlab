%%Donem Odevi -Wireless Communications ELE 561
%%Ravsan Aziz
%%Simulation of the paper "Power Allocation for Device-to-Device Communication
%%Underlaying Massive MIMO Multicasting Networks"

clear all ;  clc


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
iter_algorithm= 3 % iteration number for algorithm
Q=10; %iteration number for algorithm 2


%% Allocate positions of CUE and D2D pairs.
CUEratearray=[]; CUESNRarray=[]; D2Dratearray=[]; D2DSNRarray=[];
lambdaarray=[]; wn=[]; warray=[]; theta_array=[];

D2Dratearray2=[]; prcntg=[];
tic

D2Dratearray_varyingD=zeros(iter,10);
CUEratearray_varyingD=zeros(iter,10);
for D=1:D
    D
    for monte=1:iter
        monte
        %Ythr = 10^(Yth/10)*ones(K,1);
        rand('state',monte)
        randn('state',monte)
        
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
        wBF=sum(h,2)/sqrt(N*sum(betah.^2));
        w=wBF;
        SNRCUENoD2D=zeros(K,1);
        for k=1:K %SNR of CUEs without any D2D transmission
            SNRCUENoD2D(k)=Pbsr*abs(h(:,k)'*w)^2/No;
        end
        %pause
        for iteralgrz=1:iter_algorithm
        lambda=ones(D+1,1); %initialize power factors of D2D tx and one enodeB
        %pause
        
        while 1==1
            SNRCUE=zeros(K,1);
            InterfCUE=zeros(K,1);
            for k=1:K
                InterfCUE(k)=sum(Pdr*lambda(1:D).*(abs(g(:,k)).^2));
                SNRCUE(k)= lambda(D+1)*Pbsr*(abs(h(:,k)'*w)^2)/(InterfCUE(k)+No);
            end
            
            %         SNRCUE
            %         lambda
            %         pause
            if sum(SNRCUE<Ythr-1e-8)==0
                break;
            else
                
                [maxinterf,maxindex]=min(SNRCUE);
                Idrop=InterfCUE(maxindex)-(Pbsr*(abs(h(:,maxindex)'*w)^2)/Ythr-No);
                for d=1:D %calculate share of drop for each D2D
                    Ishare(d)=lambda(d)*Pdr*(abs(g(d,maxindex)).^2)/InterfCUE(maxindex);
                    C(d)=(abs(rho(d,d)))^2/(trace((abs(rho)).^2));
                end
                for ii=1:100
                    %ii
                    ksi=10^(-ii);
                    lambdanew=lambda;
                    for d=1:D
                        W(d)=(1-ksi)*Ishare(d)+ksi*(1-C(d));
                        lambdanew(d)=lambda(d)-W(d)*Idrop/(Pdr*(abs(g(d,maxindex)).^2)*sum(W));
                    end
                    if sum(lambdanew<0)==0
                        lambda=lambdanew;
                        break;
                    elseif sum(lambdanew<0)==length(lambda(1:D))
                        lambda(1:D)=0;
                        break;
                    end
                end
            end
        end
        
        
        for algiter=1:Q
            wfirst=wBF;
            wn=[];
            %wn(:,1)=wBF;
            wn(:,1)=w;
            %% Part 2 Optimize w
            A=zeros(1,D); B=zeros(1,D); Yk=zeros(1,K);
            for d=1:D
                A(d)=(1/(Pbsr*lambda(D+1)))*(sum(lambda(1:D).*Pdr.*(abs(rho(:,d)).^2))+No);
                B(d)=(1/(Pbsr*lambda(D+1)))*(sum(lambda(1:D).*Pdr.*(abs(rho(:,d)).^2))-lambda(d)*Pdr*(abs(rho(d,d))^2)+No);
            end
            for k=1:K
                Yk(k)=(Ythr/(Pbsr*lambda(D+1)))*(sum(lambda(1:D).*Pdr.*(abs(g(:,k)).^2))+No);
            end
            
            
            w_gradient=zeros(N,Q);
            r=zeros(N,Q);
            %tic
            for q=1:Q-1
                for d=1:D
                    w_gradient(:,q)= w_gradient(:,q)+2*(f(:,d)*f(:,d)'*wn(:,q))*(A(d)-B(d))./((abs((f(:,d))'*wn(:,q)).^2 + A(d))*(abs((f(:,d))'*wn(:,q)).^2 + B(d)));
                end
                for k=1:K
                    w_gradient(:,q)= w_gradient(:,q)-2*(h(:,k)*h(:,k)'*wn(:,q))./(abs((h(:,k))'*wn(:,q)).^2-Yk(k));
                end
                %surface projection
                r(:,q)=(w_gradient(:,q)-wn(:,q)'*w_gradient(:,q)*wn(:,q))./norm(w_gradient(:,q)-wn(:,q)'*w_gradient(:,q)*wn(:,q));
                
                %argmin
                w_theta=[];
                min_result=1e20;
                min_theta=1000;
                
                for theta=0:1:360
                    w_theta=cos(2*pi*theta/360)*wn(:,q)+sin(2*pi*theta/360)*r(:,q);
                    %f(w)
                    fw=sum(log2((abs(f'*w_theta).^2+B')./(abs(f'*w_theta).^2+A')))-sum(log2(abs(h'*w_theta).^2-Yk'));
                    %check min
                    if fw<=min_result
                        min_theta=theta;
                        min_result=fw;
                    end
                end
                
                theta_array(q)=min_theta;
                wn(:,q+1)=cos(2*pi*min_theta/360)*wn(:,q)+sin(2*pi*min_theta/360)*r(:,q);
                
                
            end
            %toc
            w=wn(:,end);
        end
            %%RATE
            SNRCUE=zeros(K,1);
            InterfCUE=zeros(K,1);
            RateCUE=zeros(K,1);
            for k=1:K %SNR of CUEs without any D2D transmission
                InterfCUE(k)=sum(Pdr*lambda(1:D).*(abs(g(:,k)).^2));
                SNRCUE(k)= lambda(D+1)*Pbsr*(abs(h(:,k)'*w)^2)/(InterfCUE(k)+No);
                RateCUE(k)=log2(1+SNRCUE(k));
            end
            %SNRCUE
            %lambda
            
%             if (sum(SNRCUE<Ythr)==0)||(lambda(1)==0)
%                 break;
%             else
%                 lambda(1:D)=max(0,lambda(1:D)-0.05);
%             end
            
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
    end
    
end


%meanpercent=mean(prcntg)
D2Dratearray_varyingD(isnan(D2Dratearray_varyingD))=0;
CUEratearray_varyingD(isnan(CUEratearray_varyingD))=0;
GD_Heuristic_D2Dsumrate_withoutprecoder=nanmean(D2Dratearray_varyingD); %average D2Dsum rate with BF precoder
CUEsumrate_withoutprecoder=nanmean(CUEratearray_varyingD); %average sum rate with BF precoder

%Plot figures
figure('Name','Sum rate of D2D pairs','NumberTitle','off')
plot(GD_Heuristic_D2Dsumrate_withoutprecoder,'-o','color','g')
grid;
xlabel('Number of D2D pairs (D)');
ylabel('Average Sum rate of D2D pairs (bits/sec)');
title('Sum rate of D2D pairs')
legend('Prec=GD PA=Uniform')

figure
cdfplot(10*log10(CUESNRarray))

toc
save(['GD_heuristic_iter' num2str(iter) '_D' num2str(D) '_' date '.mat'] )
