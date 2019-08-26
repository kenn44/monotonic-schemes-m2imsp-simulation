function explicit
%============= Variables =================
T	=10;
L	=20; %return to 200
t	=linspace(0,T,L);
dt	=t(2)-t(1);

H0	=[1 0;0 2];
H1	=[0 1;1 0];


alpha = 10;
y0 = [1;0];
ycible	=[0;1];

maxiter	= 10; %return to 30 after
c	= ones(1,L-1)*0.3; %eps
ctil = ones(1,L-1)*0.5; %eps tilde
%division by zero in mus if c=ctilde

y=zeros(2,L);
y(:,1)=y0;
p=zeros(2,L);
p(:,L)=ycible;%to edit?

psit=zeros(2,L);
psit(:,1)=expm(-(H0*dt)/(2*i))*y(:,1); %initial value

psib=zeros(2,L);
psib(:,1)=expm((H0*dt)/(2*i))*y(:,1);

chit=zeros(2,L);
chit(:,L)=expm((H0*dt)/(2*i))*p(:,L); %end value

chib=zeros(2,L);
chib(:,L)=expm(-(H0*dt)/(2*i))*p(:,L);

Jtab=zeros(1,maxiter);


%============= Main =================
for iter=1:maxiter
  %iter
  %use Newton or secante?
  %step1
    for j=L-1:-1:1
    chib(:,j+1)=expm(-(H0*dt)/(2*i))*p(:,j+1);
    etha(j)=alpha/(alpha+dt*real(chib(:,j+1).'*(H1^2)*psit(:,j+1)));
    ctil(j)=(1-etha(j))*c(j)-(etha(j)/alpha)*imag(chib(:,j+1).'*H1*psit(:,j+1));
    chit(:,j)=expm(i*H1*ctil(:,j)*dt)*chib(:,j+1);
    p(:,j)=expm(i*H0*dt/2)*chit(:,j);
  end
  %ctil
  
  
  %step2
  for j=1:L-1
    psib(:,j)=expm((H0*dt)/(2*i))*y(:,j);
    delta(j)=alpha/(alpha+dt*real(chit(:,j).'*(H1^2)*psib(:,j)));
    c(j)=(1-delta(j))*ctil(j)-(delta(j)/alpha)*imag(chit(:,j).'*H1*psib(:,j));
    psit(:,j+1)=expm(-i*H1*c(:,j)*dt)*psib(:,j);
    y(:,j+1)=expm(-i*H0*dt/2)*psit(:,j+1);
  end
  %c

  
  %J=2-norm(ycible-y(:,L))^2 - alpha * dt * norm(c)^2;
  J=2*real(ycible.'*y(:,L)) - alpha * dt * norm(c)^2;
  
  Jtab(iter)=J;
	plot(Jtab);
  xlabel("Number of iteration");
  ylabel ('{\it J_{\Delta T}(\epsilon)}')
  legend ("Explicit scheme");
	pause(.1);
  
	fprintf(2,'Iter=%i|J=%f \n',iter,J)
end
fprintf(1,'y(T)=%f \n',y(:,L))

end

%============= Functions =================
%********** Mu star **********
function mus=calculMus(x,y,mu,dt)
  Id2 = eye (size (mu));
  mus=i*(expm(mu*(x-y)*dt/i)-Id2)/(dt*(x-y));
end
