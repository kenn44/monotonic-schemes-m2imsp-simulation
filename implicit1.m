function implicit1
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

maxiter	= 30;
c	= ones(1,L-1); %epsilon
ctil = ones(1,L-1)*2; %epsilon tilde
%*2 car probleme de division par 0 dans mus si c=ctilde

Jtab=zeros(1,maxiter);

psit=zeros(2,L);
psib=zeros(2,L);
chit=zeros(2,L);
chib=zeros(2,L);

aj = ones(1,L-1);

%sinon erreur p et y undefined au debut dans psib et chit
y=zeros(2,L);
y(:,1)=y0;
p=zeros(2,L);
p(:,L)=y(:,L)-ycible;

%============= Main =================
for iter=1:maxiter
  %schema implicite, utiliser la methode de la secante?
  %step1
  
  psib = calculPsib(y,L,H0,dt); %k+1
  %compute \epsilon_{j}^{k+1}
  for j=1:L-1
    %mus1=calculMus(ctil(j),c(j),H1,dt);
    %aj(:,j) = - (1/alpha)*imag(chit(:,j)'*calculMus(ctil(j),c(j),H1,dt)*psib(:,j));
    f = @(x)x+(1/alpha)*imag(chit(:,j)'*calculMus(ctil(j),x,H1,dt)*psib(:,j));
    if (j==1)
      x0=1;
    else
      x0=c(j-1);
    endif
    x=fzero(f,x0)
    %fzero permet de chercher la racine d'une equation non lineaire
    %x represente c(j)
    c(j)=x;
  end
  chit = calculChit(p,L,H0,dt); %k
  
  %compute \psi_{j+1}^{k+1} from \psi_{j}^{k+1}
  y = calculY(y0,L,H0,H1,dt,c); %k+1
  
  
  
  %c = calculAj(ctil,c,chit,psib,H1,alpha,L,dt); 
  
  
  %step2
  %ctil=calculBj(c,ctil,chib,psit,H1,alpha,L,dt);
  %p=calculP(L,ycible,y,H0,H1,c,dt);
  %psit=calculPsit(y,L,H0,dt);
  %chib=calculChib(p,L,H0,dt);
  
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
J=2-norm(y(:,L)-ycible)^2 - alpha * dt * norm(c).^2;

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
function p=calculP(L,ycible,y,H0,H1,c,dt)
p=zeros(2,L);
p(:,L)=y(:,L)-ycible; %to edit
for j=L:-1:2
	p(:,j-1)=expm(i*H1*c(1,j-1)*dt)*expm(i*H0*dt)*p(:,j); %is it true?
end	
end
%********** Psi Tilde **********
function psit=calculPsit(y,L,H0,dt)
psit=zeros(2,L);
for j=1:L
	psit(:,j)=expm(-(H0*dt)/(2*i))*y(:,j);
end
end
%********** Psi Breve **********
function psib=calculPsib(y,L,H0,dt)
psib=zeros(2,L);
for j=1:L
	psib(:,j)=expm((H0*dt)/(2*i))*y(:,j);
end
end
%*********** Chi Tilde *********
function chit=calculChit(p,L,H0,dt)
chit=zeros(2,L);
for j=1:L
	chit(:,j)=expm((H0*dt)/(2*i))*p(:,j);
end
end
%************ Chi Breve ********
function chib=calculChib(p,L,H0,dt)
chib=zeros(2,L);
for j=1:L
	chib(:,j)=expm(-(H0*dt)/(2*i))*p(:,j);
end
end
%********** Mu star **********
function mus=calculMus(x,y,mu,dt)
  mus=i*(expm(mu*(x-y)*dt/i)-1)/(dt*(x-y));
end
%********** aj ****************
function aj=calculAj(x,y,chit,psib,mu,alpha,L,dt)
  aj = ones(1,L-1);
  mus=calculMus(x,y,mu,dt);
  for j=1:L
    aj(:,j) = - (1/alpha)*imag(chit(:,j)'*mus*psib(:,j));
  end
end
%********** bj ****************
function bj=calculBj(x,y,chib,psit,mu,alpha,L,dt)
  bj = ones(1,L-1);
  mus=calculMus(x,y,mu,dt);
  for j=1:L
    bj(:,j) = - (1/alpha)*imag(chib(:,j)'*mus*psit(:,j));
  end
end

