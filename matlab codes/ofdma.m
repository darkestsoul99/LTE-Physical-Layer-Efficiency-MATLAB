%% OFDMA 
 % Subcarrier sayisi 
 NS=512;
 % Giris sinyali olusturulmasi
 x=rand(1,NS)>0.5;
 fftlength=512;
 BW=5e6; % Bandwidth 
 FS=2*BW;%Sampling Frekansi
 pathDelays = [0 2e-5]; %comm.RayleighChannel fonksiyonunu kullanabilmek için
 pathPower = [0 -9]; %dB cinsinden. Rayleigh channel için.
 fD = 100; %Hz cinsinden. Maximum doppler shift.
 fs = 1000; %Hz cinsinden. Input sample rate.
 %gelen data'nin seriden paralele cevrimi 
 p=series2parallel(x,NS);
 % M-ary PSK & QAM Modulasyon 
 M=2;
 X=0;
 for count1=2:1:7
    if(M==2||M==4||M==16||M==64)
      M=M+X;
        %M-ary modulasyon 
      if(M<=8)
        %PSK modulasyon 
        y = pskmod(p,M);
      else
        %QAM modulasyon
        y = qammod(p,M);
      end
      ylen=length(y);
      %Mapping fonksiyonu 
      q_out=ofdma_mapping(y,ylen);
      %Inverse Fast Fourier Transform 
      outifft=ifft(q_out);
      %Cyclic Prefix eklenmesi 
      cp(count1,:)=cyclicpad(outifft,64);
      %CP uzunlugu
      cplength=length(cp);
      %Datanin paralelden seriye cevrimi 
      out=reshape(cp(count1,:),1,cplength);
      %AWGN kanali 
      ynoisy=awgn(out,100,'measured');
      % Rayleigh fading eklenmesi 
      c=comm.MIMOChannel('SampleRate',fs,'PathDelays',pathDelays, ...
          'AveragePathGains',pathPower,'MaximumDopplerShift',fD,'NumTransmitAntennas',576,'SpatialCorrelationSpecification','None');
      rf=c(ynoisy); 
      %Data'nin seriden paralele cevrimi 
      p2=series2parallel(ynoisy,cplength);
      re_par=real(p2);
      %Cyclic prefix kaldirilmasi
      rcp(count1,:)=decyclicpad(p2,64);
      rcplength=length(rcp);
      %Fast Fourier Transform
      zzfft=fft(rcp(count1,:),fftlength);
      %Demapping fonksiyonu 
      qq_out=ofdma_demapping(zzfft);
      outfft=qq_out;
      if(M<=8)
        %Receiver'da PSK Demodulasyonu 
        z=pskdemod(outfft,M);
      else
        %Receiver'da QAM Demodulasyonu 
        z=qamdemod(outfft,M);
      end
      %Data'nin paralelden seriye cevrimi
      xdash=reshape(z,1,NS);
      berr=0;
      for a=1:1:NS
        if(xdash(:,a)==x(:,a))
            berr=0;
        else
            berr=berr+1;
        end
      end
      tberr(count1,:)=berr;
      Eb_No=0:1:NS-1;
      Eb_No=0.4*Eb_No;
      if(M<=8)

ber(count1,:)=berawgn(Eb_No,'psk',M,'nondiff');
        Pe(count1,:)=erfc(sqrt(0.9*Eb_No)*sin(pi/M));
      else
        ber1(count1,:)=berawgn(0.9*Eb_No,'qam',M);
        Pe(count1,:)=2*((1-(1/sqrt(M)))*erfc(sqrt((1.5*Eb_No)/(M-1))));
      end
      for init=1:1:32
          switch M
          end
      end
   end
   M=2^count1;
end
 figure()
 %SNR ve BER Plot 

semilogy(Eb_No,ber(2,:),'k',Eb_No,ber(3,:),'g',Eb_No,ber1(5,:),'b',Eb_No,ber1(7,:),'r');
    axis([0 25 0.0001 1]);
    xlabel('SNR[dB]')
    ylabel('BER')
    legend('BPSK','QPSK','16-QAM','64-QAM')
    title('OFDMA')
    figure()
    %Error Probability Plot 
 semilogy(Eb_No,Pe(2,:),'k',Eb_No,Pe(3,:),'r',Eb_No,Pe(5,:),'b',Eb_No,Pe(7,:),'g');
    axis([0 50 0.0001 1]);
    xlabel('SNR[dB]')
    ylabel('Probability of Error')
    legend('BPSK','QPSK','16-QAM','64-QAM')
    title('OFDMA')
    h=spectrum.periodogram;
    figure()

HS=psd(h,outifft,'SpectrumType','twosided','NFFT',512,'FS',FS);
    plot(HS)
    xlabel('Sampling Frequency (2 * BW)in MHz')
    ylabel('Power Spectral Density [dBm/Hz]')
    title('OFDMA')
    grid off;