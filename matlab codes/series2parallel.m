function y = series2parallel(x,NS)
 L=length(x);
 q=floor(L/NS);
 newvec=zeros(NS,q);
 for i=1:q
 newvec(1:NS,i)=x((1+(i-1)*NS):i*NS);
 end
 y=newvec;
end
