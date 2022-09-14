function y=decyclicpad(X,L)
 N=length(X(:,1));
 Y=X(L+1:N,:);
 y=Y;
end
