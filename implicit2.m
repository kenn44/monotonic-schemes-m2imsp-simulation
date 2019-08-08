function implicit2
%============= Variables =================
T	=10;
L	=200;
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
%*2 car probleme de division par 0 dans mus si c=ctilde

y=zeros(2,L);
y(:,1)=y0;
p=zeros(2,L);
p(:,L)=y(:,L)-ycible;%to edit?

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
  iter
  %schema implicite, utiliser la methode de la secante?
  %step1
  for j=1:L-1
    psib(:,j)=expm((H0*dt)/(2*i))*y(:,j); %k+1 %check if j is ok
    %compute \eps_{j}^{k+1}
    %aj(:,j) = - (1/alpha)*imag(chit(:,j)'*calculMus(ctil(j),c(j),H1,dt)*psib(:,j));
    f = @(x)x+(1/alpha)*imag(chit(:,j)'*calculMus(ctil(j),x,H1,dt)*psib(:,j));
    if (j==1)
      x0=c(1);
    else
      x0=c(j-1);
    endif
    x=fzero(f,x0);
    %fzero permet de chercher la racine d'une equation non lineaire
    %x represente c(j)
    c(j)=x;
    chit(:,j)=expm((H0*dt)/(2*i))*p(:,j); %k
    
    %compute \psi_{j+1}^{k+1} from \psi_{j}^{k+1}
    y(:,j+1)=expm(-i*H0*dt)*expm(-i*H1*c(1,j)*dt)*y(:,j); %is it true?
  end
  c
  
  
  %step2
  for j=L-1:-1:1
    chib(:,j)=expm(-(H0*dt)/(2*i))*p(:,j); %k+1 %check if j is ok
    psit(:,j)=expm(-(H0*dt)/(2*i))*y(:,j); %k+1
    
    %bj(:,j) = - (1/alpha)*imag(chib(:,j)'*calculMus(c(j),ctil(j),H1,dt)*psit(:,j));
    f = @(x)x+(1/alpha)*imag(chib(:,j)'*calculMus(c(j),x,H1,dt)*psit(:,j));
    if (j==L-1)
      x0=ctil(L-1);
    else
      x0=ctil(j+1);
    endif
    x=fzero(f,x0);
    ctil(j)=x
    
    p(:,j)=expm(i*H1*ctil(1,j)*dt)*expm(i*H0*dt)*p(:,j+1); %is it true?
  end
  
  %J=fonctionelle(y0,L,H0,H1,dt,c,ycible,alpha);
  %Jtab(iter)=J;
	%plot(Jtab);
  %xlabel("Number of iteration");
  %ylabel ('{\it J_{\Delta T}(\epsilon)}')
  %legend ("Implicit scheme");
	%pause(.1);
  
	%fprintf(2,'Iter=%i|J=%f \n',iter,J)
end

end

%============= Functions =================
%quelle est la bonne fontionnelle ?
%quel est son gradient?
%ecrire le test debug dans une fonction
function J=fonctionelle(y0,L,H0,H1,dt,c,ycible,alpha)

y=calculY(y0,L,H0,H1,dt,c);
J=2-norm(y(:,L)-ycible)^2 - alpha * dt * norm(c)^2;

end
%********** Psi (state) **********
function y=calculY(y0,L,H0,H1,dt,c)
y=zeros(2,L);
y(:,1)=y0;
for j=2:L
	y(:,j)=expm(-i*H0*dt)*expm(-i*H1*c(1,j-1)*dt)*y(:,j-1); %is it true?
end
end
%********** Chi (adjoint state) **********
function p=calculP(L,ycible,y,H0,H1,ctil,dt)
p=zeros(2,L);
p(:,L)=y(:,L)-ycible; %to edit
for j=L:-1:2
	p(:,j-1)=expm(i*H1*ctil(1,j-1)*dt)*expm(i*H0*dt)*p(:,j); %is it true?
end	
end

%********** Mu star **********
function mus=calculMus(x,y,mu,dt)
  mus=i*(expm(mu*(x-y)*dt/i)-1)/(dt*(x-y));
end
