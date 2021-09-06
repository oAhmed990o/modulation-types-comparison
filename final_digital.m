%%
clc;
clear all;

%%
%SNR range 0 -> 60 db , 4 db step
SNR = [0 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60];
no_of_bits = 1e4;
vector = zeros(1,no_of_bits);
complex_i = sqrt(-1);

%generate random binary vector
for i = 1:no_of_bits
    vector(i) = mod(randi(2),2);
end

%%
ber_ook = zeros(1,length(SNR));
ber_prk = zeros(1,length(SNR));
ber_fsk = zeros(1,length(SNR));

for j = 1:length(SNR)

    %modulate signal according to type of modulation:
    % ook -> no change in bits required
    ook_modulated_signal = vector;

    % prk -> represent bit 1 by 1 and bit 0 by -1 ( formula : 2*vector_bits - 1)
    prk_modulated_signal =  2*vector -1;

    % fsk -> modulate 1st bit of the stream on a carrier and the other bit
    % on a carrier orthogonal on it  ( if bit_to_send = 0 send 1, else send i)
    fsk_modulated_signal = zeros(1,no_of_bits);
    
    for i = 1:no_of_bits
        if vector(i) == 0
            fsk_modulated_signal(i) = 1;
        end
        if vector(i) ~= 0
            fsk_modulated_signal(i) = complex_i;
        end
    end
    
    %snr in watt
    SNR_watt = 10^(SNR(j)/10);
    %power
    ook_power = mean(ook_modulated_signal.^2);
    prk_power = mean(prk_modulated_signal.^2);
    fsk_power = mean(abs(fsk_modulated_signal).^2);
    %noise
    ook_noise = (sqrt(ook_power/SNR_watt)/sqrt(2)).*(randn(1,no_of_bits)+complex_i*randn(1,no_of_bits));
    prk_noise = (sqrt(prk_power/SNR_watt)/sqrt(2)).*(randn(1,no_of_bits)+complex_i*randn(1,no_of_bits));
    fsk_noise = (sqrt(fsk_power/SNR_watt)/sqrt(2)).*(randn(1,no_of_bits)+complex_i*randn(1,no_of_bits));

    %apply noise to bits (symbols in fsk)
    % calculate signal power as it's not unity --> measured
    %   (rx_sequence = bits+noise)
    %or (rx_sequence = awgn(bits,SNR,'measured')
    ook_rx_sequence = ook_modulated_signal + ook_noise; 
    prk_rx_sequence = prk_modulated_signal + prk_noise;
    fsk_rx_sequence = fsk_modulated_signal + fsk_noise;

    %decide whether the rx_sequence is 1 or 0 (use relational operators and
    %indexing) 
    for i = 1: no_of_bits
        if real(ook_rx_sequence(i)) < 0.5
            ook_rx_sequence(i) = 0;
        elseif real(ook_rx_sequence(i)) >= 0.5
            ook_rx_sequence(i) = 1;
        end
        
        if real(prk_rx_sequence(i)) < 0
            prk_rx_sequence(i) = 0;
        elseif real(prk_rx_sequence(i)) >= 0
            prk_rx_sequence(i) = 1;
        end
        
        if real(fsk_rx_sequence(i)) < imag(fsk_rx_sequence(i))
            fsk_rx_sequence(i) = 1;
        elseif real(fsk_rx_sequence(i)) >= imag(fsk_rx_sequence(i))
            fsk_rx_sequence(i) = 0;
        end
       
    end
    
    %compare the original bits with the detected bits and calculate no of
    %errors (use xor or biterr)
    ook_error_bits = 0;
    prk_error_bits = 0;
    fsk_error_bits = 0;
    for i=1:no_of_bits
        if vector(i) ~= ook_rx_sequence(i)
            ook_error_bits = ook_error_bits + 1;
        end
        if vector(i) ~= prk_rx_sequence(i)
            prk_error_bits = prk_error_bits + 1;
        end
        if vector(i) ~= fsk_rx_sequence(i)
            fsk_error_bits = fsk_error_bits + 1;
        end
    end

    ook_ber =  ook_error_bits/no_of_bits;
    ber_ook(j) = ook_ber;
    
    prk_ber = prk_error_bits/no_of_bits;
    ber_prk(j) = prk_ber;
    
    fsk_ber = fsk_error_bits/no_of_bits;
    ber_fsk(j) = fsk_ber;

    %save probability of error of each SNR in matrix, ber
    % ber = [ber new prob. of error]
    ber = [ook_ber prk_ber fsk_ber];

end 

%%
%plot the ber curve against SNR (use semilogy)
figure(1)
semilogy(SNR, ber_ook);
title('modulations');grid on
hold on
semilogy(SNR, ber_prk);
hold on
semilogy(SNR, ber_fsk);
hold off
legend('OOK', 'PRK', 'FSK');

%%
%evaluate the same curves using (modem.pskmod, modem.pammod, ...)
ber_ook = zeros(1,length(SNR));
ber_prk = zeros(1,length(SNR));
ber_fsk = zeros(1,length(SNR));

for j = 1:length(SNR)

    %modulate signal according to type of modulation:
    ook_modulated_signal = vector;
    prk_modulated_signal =  pskmod(vector, 2);
    fsk_modulated_signal = fskmod(vector,2,0.5,2);
    
    %apply noise to bits (symbols in fsk)
    % calculate signal power as it's not unity --> measured 
    %   (rx_sequence = bits+noise)
    %or (rx_sequence = awgn(bits,SNR,'measured')
    
    %snr in watt
    SNR_watt = 10^(SNR(j)/10);
    %power
    ook_power = mean(ook_modulated_signal.^2);
    prk_power = mean(abs(prk_modulated_signal).^2);
    fsk_power = mean(abs(fsk_modulated_signal).^2);
    %noise
    ook_noise = (sqrt(ook_power/SNR_watt)/sqrt(2)).*(randn(1,no_of_bits)+complex_i*randn(1,no_of_bits));
    prk_noise = (sqrt(prk_power/SNR_watt)/sqrt(2)).*(randn(1,no_of_bits)+complex_i*randn(1,no_of_bits));
    fsk_noise = (sqrt(fsk_power/SNR_watt) /sqrt(2)) .* (randn(1,2*no_of_bits)+complex_i*randn(1,2*no_of_bits));
    
    %apply noise to bits (symbols in fsk)
    % calculate signal power as it's not unity --> measured
    %   (rx_sequence = bits+noise)
    %or (rx_sequence = awgn(bits,SNR,'measured')
    ook_rx_sequence = ook_modulated_signal + ook_noise; 
    prk_rx_sequence = prk_modulated_signal + prk_noise;
    fsk_rx_sequence = fsk_modulated_signal + fsk_noise;   
    
    
    %decide whether the rx_sequence is 1 or 0 (use relational operators and
    %indexing) 
    
    for i = 1: no_of_bits
        if real(ook_rx_sequence(i)) < 0.5
            ook_rx_sequence(i) = 0;
        elseif real(ook_rx_sequence(i)) >= 0.5
            ook_rx_sequence(i) = 1;
        end
    end
    
    prk_rx_sequence = pskdemod(prk_rx_sequence, 2);
    fsk_rx_sequence = fskdemod(fsk_rx_sequence,2,0.5,2);
    
    %compare the original bits with the detected bits and calculate no of
    %errors (use xor or biterr)
    ook_error_bits = 0;
    prk_error_bits = 0;
    fsk_error_bits = 0;
    for i=1:no_of_bits
        if vector(i) ~= ook_rx_sequence(i)
            ook_error_bits = ook_error_bits + 1;
        end
        if vector(i) ~= prk_rx_sequence(i)
            prk_error_bits = prk_error_bits + 1;
        end
        if vector(i) ~= fsk_rx_sequence(i)
            fsk_error_bits = fsk_error_bits + 1;
        end
    end

    ook_ber =  ook_error_bits/no_of_bits;
    ber_ook(j) = ook_ber;
    
    prk_ber = prk_error_bits/no_of_bits;
    ber_prk(j) = prk_ber;
    
    fsk_ber = fsk_error_bits/no_of_bits;
    ber_fsk(j) = fsk_ber;

    %save probability of error of each SNR in matrix, ber
    % ber = [ber new prob. of error]
    ber = [ook_ber prk_ber fsk_ber];
end

%%
%plot the matlab modulated signals
figure(2)
semilogy(SNR, ber_ook);
title('matlab mod');grid on
hold on
semilogy(SNR, ber_prk);
hold on
semilogy(SNR, ber_fsk);
hold off
legend('OOK', 'PRK', 'FSK');

%%
%evaluate the probability of error of the 16 qam modulation

ber_qam16 = zeros(1,length(SNR));

for j = 1:length(SNR)

    qam16_modulated_signal =  qammod(vector, 16);

    %apply noise to bits (symbols in fsk)
    % calculate signal power as it's not unity 
    %   (rx_sequence = bits+noise)
    %or (rx_sequence = awgn(bits,SNR,'measured')
    qam16_rx_sequence = awgn(qam16_modulated_signal,SNR(j),'measured');
    
    %decide whether the rx_sequence is 1 or 0
    qam16_rx_sequence = qamdemod(qam16_rx_sequence, 16);
    
    %compare the original bits with the detected bits and calculate no of
    %errors (use xor or biterr)
    qam16_error_bits = 0;
    for i=1:no_of_bits
        if vector(i) ~= qam16_rx_sequence(i)
            qam16_error_bits = qam16_error_bits + 1;
        end
    end

    qam16_ber =  qam16_error_bits/no_of_bits;
    ber_qam16(j) = qam16_ber;
    
end

%%
%plot the matlab qam16 modulated signal
figure(3)
%subplot(411)
semilogy(SNR, ber_qam16);
title('16 QAM pskmod');grid on