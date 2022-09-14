function [qout]=ofdma_mapping(qdata,fftlength)
 qout=zeros(fftlength,1);
 qout(2:27,:) = qdata(1:26,:);
 qout(487:end,:)= qdata(27:52,:);
end
