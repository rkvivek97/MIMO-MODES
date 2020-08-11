clear all;
M = 2; %2,4 - QPSK,16QAM;
N_src = 1e3;

MIMO_MODE = 1; %0- SISO
               %1- ANT SELECTION
               %2- MRC
               %3- ALAMOUTI
snr_db = [0:3:25];
Nt = 1;Nr = 4; %Nt = 1,2 Nr=1,2,4
N_iter = 200;

for snr_iter = 1:length(snr_db)
    for iter =1:N_iter
inp_bits = round(rand(1,N_src));
%MODULATION
switch(M)
    case 2
        const_table = [-1-1j,-1+1j,1+1j,1-1j];
        b = reshape(inp_bits,M,N_src/M);
        s = const_table(bi2de(b.')+1);
        s = s/sqrt(2);
    case 4
        const_table = [-3-3i -3-i -3+3i -3+i -1-3i -1-i -1+3i -1+i 3-3i 3-i 3+3i 3+i 1-3i 1-i 1+3i 1+i];
        b = reshape(inp_bits,N_src/M,M);
        s = const_table(bi2de(b)+1);
        s = s/sqrt(10);
end

snr_lin = 10^(snr_db(snr_iter)/10);

if MIMO_MODE == 3
    x(1,:) = reshape([s(1:2:end);-conj(s(2:2:end))],1,length(s));
    x(2,:) = reshape([s(2:2:end);conj(s(1:2:end))],1,length(s));
    x = x/sqrt(Nt);
else
    x = reshape(s,Nt,N_src/M/Nt);
end

y = zeros(Nr,length(x));

for r = 1:Nr
    for t = 1:Nt
        h(:,r,t) = sqrt(1/2)*(randn(1,1)+1j*randn(1,1));
        N0 = 1/snr_lin;
        noise = sqrt(N0/2)*(randn(1,length(x))+1j*randn(1,length(x)));

        y(r,:) = y(r,:) + h(:,r,t).'.*x(t,:) + noise;
    end
end

x_hat = [];
switch MIMO_MODE
    case 0 %SISO
        x_hat = y./h(:,1,1).';
        
    case 1  %SIMO ANTENNA SELECTION
        for i = 1:N_src/M/Nt
           [~,id] =  max(abs(h(1,:,1)));
           x_hat(i) = y(id,i)./h(1,id,1);
        end
    case 2  %MRC
        for i = 1:N_src/M/Nt
           %[~,id] =  max(abs(h(i,:,1)));
           x_hat(i) = (y(:,i).'*h(1,:,1)')/(h(i,:,1)*h(i,:,1)');
        end        
    case 3  %ALAMOUTI 
        for i=1:2:N_src/M
            h1 = [h(1,:,1).' h(1,:,2).'];
            h2 = [h(1,:,2)' -h(1,:,1)'];
            H = [h1;h2];
            scale = (sum(abs(H(:,1).^2)));
            tmp = (H')*[y(:,i);conj(y(:,i+1))]/scale;
            x_hat = [x_hat tmp(:).'];
        end
end
    
%DEMOD
op_bits = [];
for i =1:N_src/M
    [~,id] = min(x_hat(i)*sqrt(2) - const_table);
    op_bits = [op_bits de2bi(id-1,(M))];
    s_dec(i) = const_table(id);
end
s1(iter)= sum(abs(s*sqrt(2) - s_dec)>0.001)/length(s);
b1(iter) = sum(abs(op_bits-inp_bits))/N_src;
    end
    ser(snr_iter) = mean(s1);
    ber(snr_iter) = mean(b1);
end

figure;
semilogy(snr_db,ser,'-v');
title('SER vs SNR(dB) for QPSK');xlabel('SNR(dB)'),ylabel('SER');

figure;
semilogy(snr_db,ber,'-v');
title('BER vs SNR(dB) for QPSK');xlabel('SNR(dB)'),ylabel('BER');