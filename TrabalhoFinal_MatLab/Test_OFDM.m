% Simulates OFDM
EN=[-5:2:36]'+0*100; en = 10 .^(EN/10) ;
N=1024;
NSlot=1000;
CHANNEL='RAYL';
L=2; % L-th order diversity
Ts=4e-6; % Block duration
Tg=0.2*Ts; % Cyclic prefix durration
f=[-N/2:N/2-1]'/Ts; % frequencies

if (CHANNEL=='AWGN')
    alfa_med=1;tau=0;NRay=1;
elseif (CHANNEL=='REAL')
%    load pdp.940
%    load PDP_ChA.dat,pdp=PDP_ChA;
    load pdp_hipc.dat,pdp=pdp_hipc;
    tau=pdp(:,3);
    alpha_med=10 .^((pdp(:,1))/20);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
    NRay=length(alpha_med);
%    NRay=10;tau=[0:1:NRay-1]'*0.1/NRay*Ts; % Tg=0.1Ts
elseif (CHANNEL=='XTAP')
    NRay=32;
    tau=[0:NRay-1]'*Ts/N;
    alpha_med=ones(NRay,1);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
end;

Eb=21.5;
sigma=sqrt(Eb/2 ./en);
NSR=1/2 ./(en);
NEN=length(EN);
NErr=zeros(NEN,1);
for nn=1:NSlot

    %rand('state',nn*1234567); randn('state',nn*1234567);
    % This means the same channel for each slot

    if (CHANNEL=='REAL')
        Hk=zeros(N,L);
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end;
    elseif (CHANNEL=='XTAP')
        Hk=zeros(N,L);
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end;
    elseif (CHANNEL=='RRND')
        Hk=zeros(N,L);
        tau=rand(NRay,1)*Tg;
            for l=1:L
                alpha=ones(NRay,1).*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2*NRay);
                for nRay=1:NRay
                    Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
                end;
            end;
    elseif (CHANNEL=='RAYL')
        Hk=(randn(N,L)+j*randn(N,L))/sqrt(2);
    elseif (CHANNEL=='AWGN')
        Hk=ones(N,L).*exp(j*2*pi*rand(N,L));
    end;
    H2k=abs(Hk).^2;
    if (L==1) sH2k=H2k; else sH2k=sum(H2k')'; end;

    %Q2 = 2*sign(randn(N,1))+2*j*sign(randn(N,1));
    %usar qam mod
    %Ak_Tx=sign(randn(N,1))+j*sign(randn(N,1)) ; % +/-1 +/-j
    M = 256;
    bits = randi([0, 1],N*log2(M),1);

    Ak_Tx = qammod(bits,M,InputType ='bit');

    an_Tx=fftshift(ifft(fftshift(Ak_Tx)));

     %conversão para polares
     envelope_tx = abs(an_Tx);   
     phase_tx = angle(an_Tx);
     maxenv=max(envelope_tx);
     %SSPA
     Ssat=15*maxenv; %exemplo de valor
     envelope_tx_mean=mean(envelope_tx);
     p = 2;% P pode ser outro P=1, 2, 5, 10 ou 100
     satlevel = envelope_tx_mean*10^(Ssat/10);  
     A = envelope_tx./(1+(envelope_tx./satlevel).^(2*p)).^(1/(2*p));
     Fi = phase_tx;      
 
     an_Tx = A.*exp(j*Fi);
     Ak_Tx=fftshift(fft(fftshift(an_Tx)));
 
     plot(Ak_Tx,'.')






    for nEN=1:NEN
        Yk=zeros(N,L);
        for l=1:L
            Yk(:,l)=Ak_Tx.*Hk(:,l)+(randn(N,1)+j*randn(N,1))*sigma(nEN);
        end;
        YIk=0;
        for l=1:L
            YIk = YIk +Yk(:,l).*conj(Hk(:,l));
        end;
        YIk=YIk./sH2k;

        %Ak_Rx=sign(real(YIk))+j*sign(imag(YIk));
        Ak_Rx = qamdemod(YIk,M,OutputType = 'bit');

        %aux = sum( abs(real(Ak_Tx)-real(Ak_Rx)) + ... 
         %abs(imag(Ak_Tx)-imag(Ak_Rx)) ) / 2 ;
      
        aux = sum(abs(bits - Ak_Rx));

        NErr(nEN,1)=NErr(nEN,1)+aux;

    end;

    if (rem(nn,100)==0) nn, end;
end;

% BER in Rayleigh channel and L-branch diversity [Proakis]
aux=sqrt(en./(1+en));Pb_tr=0;
for l=0:L-1
    Pb_tr=Pb_tr+combin(L-1+l,l)*((1+aux)/2).^l;
end;
Pb_tr=Pb_tr.*((1-aux)/2).^L;

% BER in AWGN channel
%PbAWGN=q_x(sqrt(2*L*en));
Pb=NErr/NSlot/N/2;

%figure;
semilogy(EN, Pb, 'g-*', EN, Pb_tr, 'b:');
% Add a horizontal line at the desired BER with increased density
desired_BER = 1e-4;
desired_BER_line = ones(size(EN)) * desired_BER;

% Increase density by interpolating more points
desired_BER_line_dense = interp1(EN, desired_BER_line, linspace(min(EN), max(EN), 1000), 'linear');

hold on;
plot(linspace(min(EN), max(EN), 1000), desired_BER_line_dense, 'r--', 'LineWidth', 2);
hold off;
xlabel('E_b/N_0(dB)'), ylabel('BER')
axis([-5 40 1e-6 1]) % Adjust these limits based on your needs
%pause,clf;
