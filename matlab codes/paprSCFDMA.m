function paprSCFDMA()
   dataType = 'B-PSK'; % Modulation format.
   NS = 512; % Number of total subcarriers.
   Symbols = 16; % Data block size.
   Q = NS/Symbols; % Bandwidth spreading factor of SC-FDMA.
   BW = 5e6; % System bandwidth.
   Ts = 1/BW; % sampling rate.
   osf = 4; % Oversampling factor.
   Nsub = NS;
   Fsub = (0:Nsub-1)*BW/Nsub; % Subcarrier spacing of OFDMA.
   Runs = 1e3; % Number of iterations.
   papr1 = zeros(1,Runs); % Initialize the PAPR results for sc-fdma.
   papr3 = zeros(1,Runs); % Initialize the PAPR results for OFDMA
   for n = 1:Runs,
    % Generate random data.
    if strcmp(dataType , 'B-PSK') == 1
        tmp = round(rand(Symbols,2));
        tmp = tmp*2 - 1;
        data = (tmp(:,1) + 1i*tmp(:,2))/sqrt(2);
    elseif strcmp(dataType, '16QAM') == 1
        dataSet = [-3+3i -1+3i 1+3i 3+3i ...
                   -3+1i -1+1i 1+1i 3+1i ...
                   -3-1i -1-1i 1-1i 3-1i ...
                   -3-3i -1-3i 1-3i 3-3i];
        dataSet = dataSet / sqrt(mean(abs(dataSet).^2));
        tmp = ceil(rand(Symbols,1)*16);
        for k = 1:Symbols,
            if tmp(k) == 0
                tmp(k) = 1;
            end
            data(k) = dataSet(tmp(k));
        end
        data = data';
    end
    % Convert data to frequency domain.
    Z1 = fft(data);
    Z2 = fft(data);
    % Initialize the subcarriers.
    Y1 = zeros(NS,1);
    Y2 = zeros(NS,1);
    % Subcarrier mapping for SC-FDMA
    Y1(1:Symbols) = Z1;
    Y2(1:Symbols) = Z2;
    % Convert data back to time domain.
    y1 = ifft(Y1);
    y2 = ifft(Y2);
    % OFDMA modulation.
    % Time range of the OFDMA symbol.
    t = 0:Ts/osf:Nsub*Ts;
    y3 = 0;
    for k = 1:Symbols,
        y3= y3 + data(k)*exp(1i*2*pi*Fsub(k)*t);
    end
    % Calculate PAPR.
    papr3(n) = 10*log10(max(abs(y3).^2) /mean(abs(y3).^2));
    papr1(n) = 10*log10(max(abs(y1).^2) /mean(abs(y1).^2));
    papr2(n) = 10*log10(max(abs(y2).^2) /mean(abs(y2).^2));
 end
 % Plot CCDF.
 figure
 [N,Z3] = hist(papr3, 100);
 [N,Z1] = hist(papr1, 100);
 [N,Z2] = hist(papr2, 100);
 semilogy(Z1,1-cumsum(N)/max(cumsum(N)),'b')
 hold on
 semilogy(Z3,1-cumsum(N)/max(cumsum(N)),'black')
 hold off
 title ('PAPR of SC-FDMA and OFDMA for BPSK')
 xlabel ('PAPR[dB]')
 ylabel ('{PAPR(PAPR>PAPR0)}')
 grid off;
 % Save data.
 % save paprSCFDMA
 % QPSK:
 if strcmp(dataType, 'QPSK') == 1
    tmp = round(rand(Symbols,4));
    tmp = tmp*2 - 1;
    data = (tmp(:,1) + 1i*tmp(:,2))/sqrt(2);
 % 16-QAM:
 elseif strcmp(dataType , '16QAM') == 1
 dataSet = [-3+3i -1+3i 1+3i 3+3i ...
            -3+1i -1+1i 1+1i 3+1i ...
            -3-1i -1-1i 1-1i 3-1i ...
            -3-3i -1-3i 1-3i 3-3i];
 % 64-QAM:
 elseif strcmp(dataType , '64QAM') == 1
 dataSet = [-5+5i -1+5i 1+5i 5+5i ...
            -5+1i -1+1i 1+1i 5+1i ...
            -5-1i -1-1i 1-1i 5-1i ...
            -5-5i -1-5i 1-5i 5-5i];
    end
end