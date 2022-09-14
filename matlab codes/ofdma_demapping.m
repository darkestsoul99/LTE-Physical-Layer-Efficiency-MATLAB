function [qout]=ofdma_demapping(qdata)
 qdata = qdata';
 qout=zeros(length(qdata),1);
 qout(1:26,:)=qdata(2:27,:);
 qout(27:52,:)=qdata(487:end,:);
end