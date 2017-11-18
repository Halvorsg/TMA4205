function [black,red] = getRedBlack(M,N)

if mod(M,2) == 1
    black = 1:2:M*N;
    red = 2:2:M*N;
else
    X = reshape(1:M*N,M,N)';
    red = zeros(ceil(M*N/2),1);
    redCNT = 0;
    for i = 1:2:M
        red(redCNT+1:redCNT + length(X(2:2:end,i))) = X(2:2:end,i);
        redCNT = redCNT + length(X(2:2:end,i));
    end
    for i = 2:2:M
        red(redCNT+1:redCNT + length(X(1:2:end,i))) = X(1:2:end,i);
        redCNT = redCNT + length(X(1:2:end,i));
    end
    red = sort(red);
    black = ((1:2:N*M)+(2:2:M*N))'-red;
end 
end