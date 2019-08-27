function implicit
%============= Variables =================
T	=10;
L	=200; %return to 200
t	=linspace(0,T,L);
dt	=t(2)-t(1);

H0	=[1 0;0 2];
H1	=[0 1;1 0];

alpha = 10;
y0 = [1;0];
ycible	=[0;1];

maxiter	= 9; %return to 30 after (bug after 9)
c	= ones(1,L-1)*0.3; %eps
ctil = ones(1,L-1)*0.5; %eps tilde
%division by zero in mus if c=ctilde

y=zeros(2,L);
y(:,1)=y0;
p=zeros(2,L);
p(:,L)=ycible;

psit=zeros(2,L);
psit(:,1)=expm(-(H0*dt)/(2*i))*y(:,1); %initial value

psib=zeros(2,L);
psib(:,1)=expm((H0*dt)/(2*i))*y(:,1);

chit=zeros(2,L);
chit(:,L)=expm((H0*dt)/(2*i))*p(:,L); %end value

chib=zeros(2,L);
chib(:,L)=expm(-(H0*dt)/(2*i))*p(:,L);

Jtab=zeros(1,maxiter);

%a la 10eme iteration, fzero(c) ne fonctionne plus
%fzero: zero point is not bracketed
%fonctionne lorsque l'on remplace le produit hermitien x'*M*y
%par le produit scalaire x.'*M*y
%en fait, la fonctionnelle a deja converge depuis
%et est stable, c aussi
%essayer methode de Newton ou secante
%============= Main =================
for iter=1:maxiter
  %step1
    for j=L-1:-1:1
    chib(:,j+1)=expm(-(H0*dt)/(2*i))*p(:,j+1);
    f = @(x)x+(1/alpha)*imag(chib(:,j+1)'*calculMus(c(j),x,H1,dt)*psit(:,j+1));
    if (j==L-1)
      x0=ctil(L-1);
    else
      x0=ctil(j+1);
    endif
    x=fzero(f,x0);
    ctil(j)=x;
    chit(:,j)=expm(i*H1*ctil(:,j)*dt)*chib(:,j+1);
    
    p(:,j)=expm(i*H0*dt/2)*chit(:,j);
  end
  %ctil
  
  %step2
  for j=1:L-1
    psib(:,j)=expm((H0*dt)/(2*i))*y(:,j);
    f = @(x)x+(1/alpha)*imag(chit(:,j)'*calculMus(x,ctil(j),H1,dt)*psib(:,j));
    if (j==1)
      x0=c(1);
    else
      x0=c(j-1);
    endif
    x=fzero(f,x0);
    c(j)=x;

    psit(:,j+1)=expm(-i*H1*c(:,j)*dt)*psib(:,j);
    y(:,j+1)=expm(-i*H0*dt/2)*psit(:,j+1);
  end
  %c

  %J=2-norm(ycible-y(:,L))^2 - alpha * dt * norm(c)^2;
  J=2*real(ycible'*y(:,L)) - alpha * dt * norm(c)^2;
  
  Jtab(iter)=J;
	
  plot(Jtab, '--*');
  xlabel("Number of iteration");
  ylabel ('{\it J_{\Delta T}(\epsilon)}')
  legend ("Implicit scheme");
  
  %plot(t(1:end-1),c)
  %xlabel("Temps de contrôle");
  %ylabel ('Contrôle {\it \epsilon (t)}')
  %legend ("Champ de contrôle obtenu");
  
  %plot3(t,y(1,:),y(2,:))
  %xlabel("Temps");
  %ylabel ('{\it \psi (t)}')
  %legend ("Fonction d'onde");
  
	pause(.1);
  
	fprintf(2,'Iter=%i|J=%f \n',iter,J)
end
%fprintf(1,'y(T)=%f \n',y(:,L))
%cputime ()

end

%============= Functions =================
%********** Mu star **********
function mus=calculMus(x,y,mu,dt)
  Id2 = eye (size (mu));
  mus=i*(expm(mu*(x-y)*dt/i)-Id2)/(dt*(x-y));
end
