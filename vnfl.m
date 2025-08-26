clear 
clc
% Parameter value
M=20;
m=10^10;
%Velocity of the left wing of the crack
vnfl=[3.3,4.4,5.5,6.6,7.7];
%Velocity of the right wing of the crack
vnfr=6;
%The velocity of section m - b of the shaft
vne=5.8;
%Velocity of the vertical section of the wellbore
vew=5.6;
zetanfl=2;
zetanfr=2.5;
H=15;
anfl=1;
anfr=1.7;
bnfl=1;
bnfr=2;
cnfl=5;
cnfr=1;
dnfl=-1;
dnfr=1;
enfl=1;
enfr=1;
kanfl=5;
kanfr=-1;
wnfl=2;
wnfr=1;
D =5.2*10^(-4);
xnf=0.01;
z=0;
%the right wing of the crack
Lnfr=0.2;
%the left wing of the crack
Lnfl=0.25;
%Horizontal section length
L=1;
Aw=-1.9;
Bw=1.8;
An=10;
Bn=3;
eata=1;
for i=1:250
    t(i)=0.001*10^(9*i/99);
end
La = length(t);

% Code for  v(j)(equation (30))
V = zeros(1,M);
for j=1:M
   
    V(j)=0;
    for k=fix((j+1)/2):1:min(j,M/2)
                  V(j)=V(j) + k^(M/2) * factorial (2*k+1) / (factorial(M/2-k+1) * factorial(k+1) * factorial(k) * factorial(j-k+1) * factorial(2*k-j+1));% 改进

    end
    V(j)=(-1)^(M/2+j)*V(j);
end


%Inverse the pressure in Laplace space to real space
for k=1:5
   for i=1:La
      ft(i,k)=0;
      u=log(2)/t(i);
      for j=1:M
     s=j*u;
     %The characteristic root of the left wing of the first or the n-th crack,equation（27）-（28）
       lambda1nfl=(vnfl(k)+sqrt(vnfl(k)^2+4*D*(s+zetanfl)))/(2*D);
       lambda2nfl=(vnfl(k)-sqrt(vnfl(k)^2+4*D*(s+zetanfl)))/(2*D);
       %The characteristic root parameters of the left wing of the first or mth crack
       lambda1nfla=(vnfl(k)+sqrt(vnfl(k)^2+4*D*(s+zetanfl-anfl/eata)))/(2*D);
       lambda2nfla=(vnfl(k)-sqrt(vnfl(k)^2+4*D*(s+zetanfl-anfl/eata)))/(2*D);
       %The characteristic roots of the right wing of the first or the n-th crack,equation（36）-（37） 
       lambda1nfr=(vnfr+sqrt(vnfr^2+4*D*(s+zetanfr)))/(2*D);
       lambda2nfr=(vnfr-sqrt(vnfr^2+4*D*(s+zetanfr)))/(2*D);
       %Characteristic root parameters of the right-wing features of the first or n cracks-anfl/eata，
       lambda1nfra=(vnfr+sqrt(vnfr^2+4*D*(s+zetanfr-anfr/eata)))/(2*D);
       lambda2nfra=(vnfr-sqrt(vnfr^2+4*D*(s+zetanfr-anfr/eata)))/(2*D);
       %The characteristic root from the first or the n-th crack to the bottom end of the vertical wellbore,equation（55）-（56） 
       lambda1ne=(vne+sqrt(vne^2+4*D*s))/(2*D);
       lambda2ne=(vne-sqrt(vne^2+4*D*s))/(2*D);
       %The characteristic root from the bottom of the vertical wellbore to the wellhead,quation（60）-（61） 
       lambda1ew=(vew+sqrt(vew^2+4*D*s))/(2*D);
       lambda2ew=(vew-sqrt(vew^2+4*D*s))/(2*D);
       %equation（64）
       znfl=-(1-1/(2*m))*Lnfl;%(16)
       gnfl=cnfl*xnf./L+dnfl*znfl/Lnfl+wnfl;
       %equation（66）
       znfr=(1-1/(2*m))*Lnfr;%(17)
       gnfr=cnfr*xnf./L+dnfr*znfr/Lnfr+wnfr;
       %equation（69）beta2nl(s)
       beta2nl=((exp(gnfl)*((kanfl*(-lambda1nfla*znfl-1)/(lambda1nfla^2)-(bnfl*xnf+enfl)/lambda1nfla))/(lambda1nfla-lambda2nfla)...
              +exp(gnfl)*((kanfl*(-lambda2nfla*znfl-1)/(lambda2nfla^2)-(bnfl*xnf+enfl)/lambda2nfla))/(lambda2nfla-lambda1nfla))...
              +znfl*(lambda1nfl*((kanfl*(-lambda1nfl*znfl-1)/(lambda1nfl^2)-(bnfl*xnf+enfl)/lambda1nfl))/(lambda1nfl-lambda2nfl)...
              +lambda2nfl*((kanfl*(-lambda2nfl*znfl-1)/(lambda2nfl^2)-(bnfl*xnf+enfl)/lambda2nfl))/(lambda2nfl-lambda1nfl)))...
              *(-znfl*lambda2nfl*exp(lambda2nfl*znfl)+exp(lambda2nfla*znfl+gnfl-anfl/eata))^(-1);

       %equation（70）betanr2(s)
       beta2nr=((exp(gnfr)*((kanfr*(-lambda1nfra*znfr-1)/(lambda1nfra^2)-(bnfr*xnf+enfr)/lambda1nfra))/(lambda1nfra-lambda2nfra)...
              +exp(gnfr)*((kanfr*(-lambda2nfra*znfr-1)/(lambda2nfra^2)-(bnfr*xnf+enfr)/lambda2nfra))/(lambda2nfra-lambda1nfra))...
              +znfr*(lambda1nfr*((kanfr*(-lambda1nfr*znfr-1)/(lambda1nfr^2)-(bnfr*xnf+enfr)/lambda1nfr))/(lambda1nfr-lambda2nfr)...
              +lambda2nfr*((kanfr*(-lambda2nfr*znfr-1)/(lambda2nfr^2)-(bnfr*xnf+enfr)/lambda2nfr))/(lambda2nfr-lambda1nfr)))...
              *(-znfr*lambda2nfr*exp(lambda2nfr*znfr)+exp(lambda2nfra*znfr+gnfr-anfr/eata))^(-1);

       %Concentration on the left wing of the first or n-th crack,equation（67）
       Cnfl=beta2nl*exp(lambda2nfl*z)+(kanfl*(-lambda1nfl*z-1)/(lambda1nfl^2)-(bnfl*xnf+enfl)/lambda1nfl)/(lambda2nfl-lambda1nfl)...
             +(kanfl*(-lambda2nfl*z-1)/(lambda2nfl^2)-(bnfl*xnf+enfl)/lambda2nfl)/(lambda1nfl-lambda2nfl);
       %The concentration on the right wing of the first or n-th crack,equation（68）
       Cnfr=beta2nr*exp(lambda2nfr*z)+(kanfr*(-lambda1nfr*z-1)/(lambda1nfr^2)-(bnfr*xnf+enfr)/lambda1nfr)/(lambda2nfr-lambda1nfr)...
             +(kanfr*(-lambda2nfr*z-1)/(lambda2nfr^2)-(bnfr*xnf+enfr)/lambda2nfr)/(lambda1nfr-lambda2nfr);
       %Concentration at the bottom of the wellbore,equation（73）
       Cne=exp(lambda2ne*L)*((Cnfl+ Cnfr)/2+(An+Bn*lambda1ne)/((lambda2ne-lambda1ne)*lambda1ne^2)+(An+Bn*lambda2ne)/((lambda1ne-lambda2ne)*lambda2ne^2))...
           -(An*lambda1ne*L+An+Bn*lambda1ne)/((lambda2ne-lambda1ne)*lambda1ne^2)-(An*lambda2ne*L+An+Bn*lambda2ne)/((lambda1ne-lambda2ne)*lambda2ne^2);
       %Concentration of pollutant at the wellhead in Laplace-space ,equation（74）
       Cew=exp(lambda2ew*H)*(Cne+(Aw+Bw*lambda1ew)/((lambda2ew-lambda1ew)*lambda1ew^2)+(Aw+Bw*lambda2ew)/((lambda1ew-lambda2ew)*lambda2ew^2))...
           -(Aw*lambda1ew*H+Aw+Bw*lambda1ew)/((lambda2ew-lambda1ew)*lambda1ew^2)-(Aw*lambda2ew*L+Aw+Bw*lambda2ew)/((lambda1ew-lambda2ew)*lambda2ew^2);
       %Concentration of pollutant at the wellhead in real-space, equation(75)
       ft(i,k) =  ft(i,k) + V(j) *(Cew);
      end
     ft(i,k)=ft(i,k)*u;
   end
end
% Double logarithmic curves of concentration of pollutant at the wellhead in real-space at different vnfl
loglog(t,ft(:,1),t,ft(:,2),t,ft(:,3),t,ft(:,4),t,ft(:,5), 'k-','Linewidth',1.)
hold on 
axis([10^(0) 10^(4) 580 10^(4.7)]);

xlabel({'$$t (h)$$ ','','Figure 9 Concentration curves of the effect of $$v_{nfl}$$'},'FontName','Times New Roman','FontSize',10,'Interpreter','latex') 
ylabel('$$\mu_{ew} (mg/L)$$' ,'FontName','Times New Roman','FontSize',10,'Interpreter','latex')
legend('$$ v_{nfl}=3.3 $$','$$ v_{nfl}=4.4 $$','$$ v_{nfl}=5.5 $$','$$ v_{nfl}=6.6 $$','$$ v_{nfl}=7.7 $$','Interpreter','latex')
grid on  
