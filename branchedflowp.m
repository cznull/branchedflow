function [r,par]=branchedflowp(corr,n)
m=rand(n)-0.5+1i*(rand(n)-0.5);
x=[(0:n/2-1),(n/2:-1:1)]';
m=(n*n/corr)*m.*(exp(-(1.0/corr/corr)*x.*x)*exp(-(1.0/corr/corr)*x.*x)');
r=real(ifft2(m));

V=r*0.05;
u=(1:n)'-n/2-0.5;
u=exp(-(16.0/8192)*u.*u);%incident light

par=zeros(n,n);

k=2*pi*4096;%4096: 1/wave length
T=(2*pi*pi/k)*x.*x;
for i=1:n
    u=fft(u);
    u=u.*exp(1i*(1.0/n)*T);
    u=ifft(u);
    u=u.*exp(-1i*(k/n/2)*V(:,i));
    par(:,i)=u;
end
imshow(abs(par),[0,2]);
end

%example: branchedflowp(64,4096);
%4096: grid size
%64: 1/correlation length
