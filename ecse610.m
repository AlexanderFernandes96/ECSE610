%% Setup Channel
close ALL; clear ALL; 
for MQAM = [16 64 128]
    for K_subcarriers = [64 256 1024]
        fprintf('%i-QAM, %i subcarriers\n', MQAM, K_subcarriers); 
        close ALL;
        K=K_subcarriers; CP=K/8; Nofdm=K+CP;
        Nps=2; Np=K/Nps; % Pilot spacing and number of pilots per OFDM symbol 
        M=MQAM; % Number of bits per (modulated) symbol
        Es=1; Beta=sqrt(3/2/(M-1)*Es); % Signal energy and QAM normalization factor 
        L = 5; % channel tap length
        sigma = 0.5;
        mimo_width = 10; mimo_height = 10;
        Nt = 50; Nr = mimo_width * mimo_height;
        Tx = 1; Rx = 1; % Power Delay Profile Channel
        SNR = -10:5:20;

        % Rayleigh Uncorrelated Channel
        h_unc = sqrt(sigma) * (randn(L, Nr, Nt) +  1j * randn(L, Nr, Nt)) / sqrt(L);
        % Rician Channel
        K_ric = 10^( 6 / 10);
        h_ric = sqrt(K_ric/(K_ric+1)) + sqrt(1/(K_ric+1))*h_unc;
        % Correlated Channel
        fc = 3E9; % carrier frequency
        lam = 3E8 / fc; % speed of light / carrier frequency
        d = lam / 2; % distance between antennas
        h_cor = zeros(L, Nr, Nt);
        for r = 1:Nt % per UE
            Rtx_r = zeros(Nr, L);
            v_r = diag(sqrt(sigma) * (randn(L,1) +  1j * randn(L,1)));
            theta = ones(L,1)* (r/Nt) * pi; % azimuth
            phi = (-35 + 10*rand(L,1)) * pi / 180; % elevation
            for l = 1:L
                w = 0:mimo_width-1;
                h = 0:mimo_height-1;
                a_theta = exp(1j * 2*pi * w * (d/lam) * sin(theta(l)));
                a_phi = exp(1j * 2*pi * h * (d/lam) * sin(phi(l)));
                Rtx_r(:,l) = reshape(kron(a_theta.', a_phi), 1, Nr) / sqrt(L);
            end
            h_cor(:, :, r) = (Rtx_r * v_r).';
        end

        % Power Delay Profile
        Auto_h_unc = zeros(L, Nr, Nt);
        Auto_h_unc_dB = zeros(L, Nr, Nt);
        Auto_h_ric = zeros(L, Nr, Nt);
        Auto_h_ric_dB = zeros(L, Nr, Nt);
        Auto_h_cor = zeros(L, Nr, Nt);
        Auto_h_cor_dB = zeros(L, Nr, Nt);
        for q = 1:Nr
            for r = 1:Nt
                c = xcorr(h_unc(:,q,r), h_unc(:,q,r));
                c = c(L:2*L-1);
                Auto_h_unc(:, q, r) = c;
                Auto_h_unc_dB(:, q, r) = 20*log10( real(c).^2 + imag(c).^2 );

                c = xcorr(h_ric(:,q,r), h_ric(:,q,r));
                c = c(L:2*L-1);
                Auto_h_ric(:, q, r) = c;
                Auto_h_ric_dB(:, q, r) = 20*log10( real(c).^2 + imag(c).^2 );

                c = xcorr(h_cor(:,q,r), h_cor(:,q,r));
                c = c(L:2*L-1);
                Auto_h_cor(:, q, r) = c;
                Auto_h_cor_dB(:, q, r) = 20*log10( real(c).^2 + imag(c).^2 );
            end
        end

        F = dftmtx(K) ./ sqrt(K);
        f = F(:, 1:L);

        % Send Signal
        B = zeros(K, Nt);
        S = zeros(K, Nt);
        msg=randi([0, M-1], K - Np, Nt); % bit generation 
        msg_bin = de2bi(msg,'left-msb');
        ip = 0;
        for k=1:K 
            if mod(k,Nps)==1 
                % Pilot sequence generation
                B(k, :) = 2*(randn(1, Nt)>0)-1; % Random
                ip = ip + 1;
            else
                S(k, :) = Beta*qammod(msg(k - ip, :), M, 'gray'); % Add message
            end
        end
        X = S + B;
        x = ifft(X,K);
        xt = [x(K-CP+1:K,:); x]; % Add CP 

        figure(1)
        plot(0:L-1, Auto_h_unc_dB(:,Rx,Tx), '--x', 'DisplayName', 'Rayleigh')
        hold on;
        plot(0:L-1, Auto_h_ric_dB(:,Rx,Tx), '--o', 'DisplayName', 'Rician')
        plot(0:L-1, Auto_h_cor_dB(:,Rx,Tx), '--^', 'DisplayName', 'Correlation')
        hold off;
        title('Channel Power Delay Profile across channel Tx ' + ...
            string(Tx) + ' to Rx ' + string(Rx));
        ylabel('Power [dB]')
        xlabel('Delay [tap]')
        grid on;
        legend

        lag = 0;

        figure(2)
        heatmap(reshape(abs(sum(h_unc, 1)), [Nr,Nt]).')
        title('Amplitude Matrix Rayleigh Channel');
        xlabel('Rx'); ylabel('Tx');
        figure(3)
        heatmap(reshape(abs(sum(h_ric, 1)), [Nr,Nt]).')
        title('Amplitude Matrix Rician Channel');
        xlabel('Rx'); ylabel('Tx');
        figure(4)
        heatmap(reshape(abs(sum(h_cor, 1)), [Nr,Nt]).')
        title('Amplitude Matrix Correlated Channel (width\timesheight) = (' + ...
            string(mimo_width) + '\times' + string(mimo_height) + ')');
        xlabel('Rx'); ylabel('Tx');

        %% Channel Estimation
        % Choose h
        h_name = "Rayleigh"; h = h_unc; h_mrk = '--x';
        % h_name = "Rician"; h = h_ric; h_mrk = '--o';
        % h_name = "Correlation " + string(mimo_width) + "\times" + string(mimo_height); h = h_cor; h_mrk = '--^';
        H = fft(h, K);
        BER_MF = zeros(1,length(SNR));
        BER_ZF = zeros(1,length(SNR));
        BER_MMSE = zeros(1,length(SNR));
        MSE = zeros(1,length(SNR));
        for snr = 1:length(SNR)
            fprintf('SNR %3i, Rx ', SNR(snr));
            % Transmit over channel
            Y = zeros(K, Nr);
            H_est = zeros(K, Nr, Nt);
            h_est = zeros(L, Nr, Nt);
            for q = 1:Nr
                if mod(q, floor(Nr/10)) == 1
                    fprintf(' %3i',q);
                end
                y_channel = conv(xt(:,1), h(:,q,1)); % Channel path (convolution) 
                for r = 2:Nt
                    y_channel = y_channel + conv(xt(:,r), h(:,q,r));
                end

                yt = awgn(y_channel, SNR(snr), 'measured'); 
                y = yt(CP+1:Nofdm,:); Y(:,q) = fft(y, K); % Remove CP and FFT 

                % Least Squares Channel Estimation
                A = zeros(K, L*Nt);
                for r = 1:Nt
                    A(:, (r-1)*L+1:r*L) = diag(B(:, r)) * f .* sqrt(K);
                end

                pA_Y = pinv(A) * Y(:,q);
                h_q_est = zeros(L, Nt);
                for r = 1:Nt
                    h_q_est(:,r) = pA_Y((r-1)*L+1:r*L);
                end

                H_est(:,q,:) = fft(h_q_est, K);
                h_est(:,q,:) = h_q_est;
            end
            fprintf(' Estimation,');

            % Operate on each orthogonal subcarrier
            ip = 0;
            Data_extracted_MF = zeros(K-Np, Nt);
            Data_extracted_ZF = zeros(K-Np, Nt);
            Data_extracted_MMSE = zeros(K-Np, Nt);

            for k=1:K 
                % MSE 
                H_k = zeros(Nr, Nt); H_est_k = zeros(Nr, Nt);
                for r = 1:Nt
                    for q = 1:Nr
                        H_k(q,r) = H(k,q,r);
                        H_est_k(q,r) = H_est(k,q,r);
                    end
                end
                MSE(snr) = MSE(snr) + norm(H_k - H_est_k)^2 / (K*Nr*Nt);


                if mod(k,Nps)==1
                    ip=ip+1; 
                else
                    % Matched Filter Detector
                    Data_extracted_MF(k-ip, :) = H_est_k' * Y(k,:).'; 

                    % Zero Forcing Equalizer Detector
                    Data_extracted_ZF(k-ip, :) = (H_est_k'*H_est_k) \ H_est_k' * Y(k,:).'; 

                    % Minimum Mean Square Error Detector
                    var_n = 1 / 10^(SNR(snr)*0.1);
                    T_MMSE = (H_est_k'*H_est_k + 2*var_n*eye(Nt)) \ H_est_k'; 

                    Data_extracted_MMSE(k-ip, :) = T_MMSE * Y(k,:).';
                end

            end
            fprintf(' Detection\n');

            % MF Detector
            msg_detected_MF = qamdemod(Data_extracted_MF/Beta, M, 'gray'); 
            msg_bin_est_MF = de2bi(msg_detected_MF,'left-msb');
            be_MF = sum(msg_bin_est_MF(:)~=msg_bin(:));
            BER_MF(snr) = be_MF / length(msg_bin(:));

            % ZF Detector
            msg_detected_ZF = qamdemod(Data_extracted_ZF/Beta, M, 'gray'); 
            msg_bin_est_ZF = de2bi(msg_detected_ZF,'left-msb');
            be_ZF = sum(msg_bin_est_ZF(:)~=msg_bin(:));
            BER_ZF(snr) = be_ZF / length(msg_bin(:));

            % MMSE Detector
            msg_detected_MMSE = qamdemod(Data_extracted_MMSE/Beta, M, 'gray'); 
            msg_bin_est_MMSE = de2bi(msg_detected_MMSE,'left-msb');
            be_MMSE = sum(msg_bin_est_MMSE(:)~=msg_bin(:));
            BER_MMSE(snr) = be_MMSE / length(msg_bin(:));
        end

        figure(5)
        semilogy(SNR, MSE, 'DisplayName', h_name);
        hold on;
        title('Mean Square Error of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('||H_{est} - H||^2')
        xlabel('SNR [dB]')
        grid on;
        legend

        figure(6)
        semilogy(SNR, BER_MF, h_mrk, 'DisplayName', 'MF Detector ' + h_name);
        hold on
        title(h_name + 'Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend
        figure(7)
        semilogy(SNR, BER_ZF, h_mrk, 'DisplayName', 'ZF Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend
        figure(8)
        semilogy(SNR, BER_MMSE, h_mrk, 'DisplayName', 'MMSE Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend

        %% Channel Estimation
        % Choose h
        % h_name = "Rayleigh"; h = h_unc; h_mrk = '--x';
        h_name = "Rician"; h = h_ric; h_mrk = '--o';
        % h_name = "Correlation " + string(mimo_width) + "\times" + string(mimo_height); h = h_cor; h_mrk = '--^';
        H = fft(h, K);
        BER_MF = zeros(1,length(SNR));
        BER_ZF = zeros(1,length(SNR));
        BER_MMSE = zeros(1,length(SNR));
        MSE = zeros(1,length(SNR));
        for snr = 1:length(SNR)
            fprintf('SNR %3i, Rx ', SNR(snr));
            % Transmit over channel
            Y = zeros(K, Nr);
            H_est = zeros(K, Nr, Nt);
            h_est = zeros(L, Nr, Nt);
            for q = 1:Nr
                if mod(q, floor(Nr/10)) == 1
                    fprintf(' %3i',q);
                end
                y_channel = conv(xt(:,1), h(:,q,1)); % Channel path (convolution) 
                for r = 2:Nt
                    y_channel = y_channel + conv(xt(:,r), h(:,q,r));
                end

                yt = awgn(y_channel, SNR(snr), 'measured'); 
                y = yt(CP+1:Nofdm,:); Y(:,q) = fft(y, K); % Remove CP and FFT 

                % Least Squares Channel Estimation
                A = zeros(K, L*Nt);
                for r = 1:Nt
                    A(:, (r-1)*L+1:r*L) = diag(B(:, r)) * f .* sqrt(K);
                end

                pA_Y = pinv(A) * Y(:,q);
                h_q_est = zeros(L, Nt);
                for r = 1:Nt
                    h_q_est(:,r) = pA_Y((r-1)*L+1:r*L);
                end

                H_est(:,q,:) = fft(h_q_est, K);
                h_est(:,q,:) = h_q_est;
            end
            fprintf(' Estimation,');

            % Operate on each orthogonal subcarrier
            ip = 0;
            Data_extracted_MF = zeros(K-Np, Nt);
            Data_extracted_ZF = zeros(K-Np, Nt);
            Data_extracted_MMSE = zeros(K-Np, Nt);

            for k=1:K 
                % MSE 
                H_k = zeros(Nr, Nt); H_est_k = zeros(Nr, Nt);
                for r = 1:Nt
                    for q = 1:Nr
                        H_k(q,r) = H(k,q,r);
                        H_est_k(q,r) = H_est(k,q,r);
                    end
                end
                MSE(snr) = MSE(snr) + norm(H_k - H_est_k)^2 / (K*Nr*Nt);


                if mod(k,Nps)==1
                    ip=ip+1; 
                else
                    % Matched Filter Detector
                    Data_extracted_MF(k-ip, :) = H_est_k' * Y(k,:).'; 

                    % Zero Forcing Equalizer Detector
                    Data_extracted_ZF(k-ip, :) = (H_est_k'*H_est_k) \ H_est_k' * Y(k,:).'; 

                    % Minimum Mean Square Error Detector
                    var_n = 1 / 10^(SNR(snr)*0.1);
                    T_MMSE = (H_est_k'*H_est_k + 2*var_n*eye(Nt)) \ H_est_k'; 

                    Data_extracted_MMSE(k-ip, :) = T_MMSE * Y(k,:).';
                end

            end
            fprintf(' Detection\n');

            % MF Detector
            msg_detected_MF = qamdemod(Data_extracted_MF/Beta, M, 'gray'); 
            msg_bin_est_MF = de2bi(msg_detected_MF,'left-msb');
            be_MF = sum(msg_bin_est_MF(:)~=msg_bin(:));
            BER_MF(snr) = be_MF / length(msg_bin(:));

            % ZF Detector
            msg_detected_ZF = qamdemod(Data_extracted_ZF/Beta, M, 'gray'); 
            msg_bin_est_ZF = de2bi(msg_detected_ZF,'left-msb');
            be_ZF = sum(msg_bin_est_ZF(:)~=msg_bin(:));
            BER_ZF(snr) = be_ZF / length(msg_bin(:));

            % MMSE Detector
            msg_detected_MMSE = qamdemod(Data_extracted_MMSE/Beta, M, 'gray'); 
            msg_bin_est_MMSE = de2bi(msg_detected_MMSE,'left-msb');
            be_MMSE = sum(msg_bin_est_MMSE(:)~=msg_bin(:));
            BER_MMSE(snr) = be_MMSE / length(msg_bin(:));
        end

        figure(5)
        semilogy(SNR, MSE, 'DisplayName', h_name);
        hold on;
        title('Mean Square Error of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('||H_{est} - H||^2')
        xlabel('SNR [dB]')
        grid on;
        legend

        figure(6)
        semilogy(SNR, BER_MF, h_mrk, 'DisplayName', 'MF Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend
        figure(7)
        semilogy(SNR, BER_ZF, h_mrk, 'DisplayName', 'ZF Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend
        figure(8)
        semilogy(SNR, BER_MMSE, h_mrk, 'DisplayName', 'MMSE Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend

        %% Channel Estimation
        % Choose h
        % h_name = "Rayleigh"; h = h_unc; h_mrk = '--x';
        % h_name = "Rician"; h = h_ric; h_mrk = '--o';
        h_name = "Correlation " + string(mimo_width) + "\times" + string(mimo_height); h = h_cor; h_mrk = '--^';
        H = fft(h, K);
        BER_MF = zeros(1,length(SNR));
        BER_ZF = zeros(1,length(SNR));
        BER_MMSE = zeros(1,length(SNR));
        MSE = zeros(1,length(SNR));
        for snr = 1:length(SNR)
            fprintf('SNR %3i, Rx ', SNR(snr));
            % Transmit over channel
            Y = zeros(K, Nr);
            H_est = zeros(K, Nr, Nt);
            h_est = zeros(L, Nr, Nt);
            for q = 1:Nr
                if mod(q, floor(Nr/10)) == 1
                    fprintf(' %3i',q);
                end
                y_channel = conv(xt(:,1), h(:,q,1)); % Channel path (convolution) 
                for r = 2:Nt
                    y_channel = y_channel + conv(xt(:,r), h(:,q,r));
                end

                yt = awgn(y_channel, SNR(snr), 'measured'); 
                y = yt(CP+1:Nofdm,:); Y(:,q) = fft(y, K); % Remove CP and FFT 

                % Least Squares Channel Estimation
                A = zeros(K, L*Nt);
                for r = 1:Nt
                    A(:, (r-1)*L+1:r*L) = diag(B(:, r)) * f .* sqrt(K);
                end

                pA_Y = pinv(A) * Y(:,q);
                h_q_est = zeros(L, Nt);
                for r = 1:Nt
                    h_q_est(:,r) = pA_Y((r-1)*L+1:r*L);
                end

                H_est(:,q,:) = fft(h_q_est, K);
                h_est(:,q,:) = h_q_est;
            end
            fprintf(' Estimation,');

            % Operate on each orthogonal subcarrier
            ip = 0;
            Data_extracted_MF = zeros(K-Np, Nt);
            Data_extracted_ZF = zeros(K-Np, Nt);
            Data_extracted_MMSE = zeros(K-Np, Nt);

            for k=1:K 
                % MSE 
                H_k = zeros(Nr, Nt); H_est_k = zeros(Nr, Nt);
                for r = 1:Nt
                    for q = 1:Nr
                        H_k(q,r) = H(k,q,r);
                        H_est_k(q,r) = H_est(k,q,r);
                    end
                end
                MSE(snr) = MSE(snr) + norm(H_k - H_est_k)^2 / (K*Nr*Nt);


                if mod(k,Nps)==1
                    ip=ip+1; 
                else
                    % Matched Filter Detector
                    Data_extracted_MF(k-ip, :) = H_est_k' * Y(k,:).'; 

                    % Zero Forcing Equalizer Detector
                    Data_extracted_ZF(k-ip, :) = (H_est_k'*H_est_k) \ H_est_k' * Y(k,:).'; 

                    % Minimum Mean Square Error Detector
                    var_n = 1 / 10^(SNR(snr)*0.1);
                    T_MMSE = (H_est_k'*H_est_k + 2*var_n*eye(Nt)) \ H_est_k'; 

                    Data_extracted_MMSE(k-ip, :) = T_MMSE * Y(k,:).';
                end

            end
            fprintf(' Detection\n');

            % MF Detector
            msg_detected_MF = qamdemod(Data_extracted_MF/Beta, M, 'gray'); 
            msg_bin_est_MF = de2bi(msg_detected_MF,'left-msb');
            be_MF = sum(msg_bin_est_MF(:)~=msg_bin(:));
            BER_MF(snr) = be_MF / length(msg_bin(:));

            % ZF Detector
            msg_detected_ZF = qamdemod(Data_extracted_ZF/Beta, M, 'gray'); 
            msg_bin_est_ZF = de2bi(msg_detected_ZF,'left-msb');
            be_ZF = sum(msg_bin_est_ZF(:)~=msg_bin(:));
            BER_ZF(snr) = be_ZF / length(msg_bin(:));

            % MMSE Detector
            msg_detected_MMSE = qamdemod(Data_extracted_MMSE/Beta, M, 'gray'); 
            msg_bin_est_MMSE = de2bi(msg_detected_MMSE,'left-msb');
            be_MMSE = sum(msg_bin_est_MMSE(:)~=msg_bin(:));
            BER_MMSE(snr) = be_MMSE / length(msg_bin(:));
        end

        figure(5)
        semilogy(SNR, MSE, 'DisplayName', h_name);
        hold on;
        title('Mean Square Error of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('||H_{est} - H||^2')
        xlabel('SNR [dB]')
        grid on;
        legend

        figure(6)
        semilogy(SNR, BER_MF, h_mrk, 'DisplayName', 'MF Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend
        figure(7)
        semilogy(SNR, BER_ZF, h_mrk, 'DisplayName', 'ZF Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend
        figure(8)
        semilogy(SNR, BER_MMSE, h_mrk, 'DisplayName', 'MMSE Detector ' + h_name);
        hold on
        title('Bit Error Rate of LSCE Nr=' + string(Nr) + ' Nt=' + string(Nt))
        ylabel('BER')
        xlabel('SNR [dB]')
        grid on;
        legend

        subFolderName = "\" + string(K) + "K_" + string(M) + "QAM_" + string(Nt) + "Nt_" + string(mimo_width) + "w" + string(mimo_height) + "hNr\";
        if ~exist(pwd + subFolderName, 'dir')
               mkdir(pwd + subFolderName)
        end
        fprintf("saving in: %s\n", subFolderName);
        saveas(figure(1), pwd + subFolderName + "PDP.fig");
        saveas(figure(2), pwd + subFolderName + "AM_ray.fig");
        saveas(figure(3), pwd + subFolderName + "AM_ric.fig");
        saveas(figure(4), pwd + subFolderName + "AM_cor.fig");
        saveas(figure(5), pwd + subFolderName + "MSE.fig");
        saveas(figure(6), pwd + subFolderName + "BER_MF.fig");
        saveas(figure(7), pwd + subFolderName + "BER_ZF.fig");
        saveas(figure(8), pwd + subFolderName + "BER_MMSE.fig");
    end
end