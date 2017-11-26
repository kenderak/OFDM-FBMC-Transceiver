function [ H, G ] = separate_single( Y )

[N, numOfSym] = size(Y);

R = real(Y);
I = imag(Y);

H = single(complex(zeros(N,numOfSym)));
G = single(complex(zeros(N,numOfSym)));

for symId=1:numOfSym
    for n=1:N
        H(n,symId) = (R(n,symId)/2)+R(mod(N+1-n,N)+1,symId)/2+1i*((I(n,symId)/2)-(I(mod(N+1-n,N)+1,symId)/2));

        G(n,symId) = (I(n,symId)/2)+I(mod(N+1-n,N)+1,symId)/2-1i*((R(n,symId)/2)-(R(mod(N+1-n,N)+1,symId)/2));
    end
end

%H = H/N;
%G = G/N;
end


