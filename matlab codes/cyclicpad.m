function y=cyclicpad(X,L)
 N=length(X(:,1));
 %N-L+1
 Y=[X(N-L+1:N,:);X];
 y=Y;
end