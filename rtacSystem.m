clear all;
%nilai paramter 
clear all;
 M=1.3608;
 m=0.096;
 R=0.0592;
 I=2.175*10^-4;
 k=186.3;
 N=0;
 F=30;
 
 %pembagi di persamaan
 h=I*(M+m)+M*m*R^2;
 
 %nilai di matriks A
 a21=-1*k*(I+m*R^2)/h;
 a23=(((m*R*N*h)-(2*m^3*R^3*N))-(2*F*h*(I+m*R^2)*(m^2*R^2)))/(h^2);
 a41=m*R*k/h;
 a43=(-m*R*(-F+2*m*R*(M+m)*N))/(h^2);
 
 %nilai di matriks B
 b21=(I+m*R^2)/h;
 b22=(-m*R/h);
 b41=(-m*R/h);
 b42=(M+m)/h;
 
 %matriks A,B,C, dan D
 A=[0 1 0 0;a21 0 a23 0; 0 0 0 1;a41 0 a43 0];
 B=[0 0;b21 b22;0 0;b41 b42];
 C=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
 D=zeros(4,2);
 
 [n,m]=size(A);
 %================================================================
 %observable
 Ob=obsv(A,C);
 unob = n-rank(Ob);
 if(unob==0)
    disp('Given System is Observable.');
 else
    disp('Given System is Unobservable');
 end
 
 %controllable
 Co = ctrb(A,B);
 unco=m-rank(Co);
 if(unco==0)
    disp('Given System is Controllable.');
 else
    disp('Given System is Uncontrollable');
 end
 
 %================================================
 %Mengecek kestabilan (melihat apakah pole berada di half-right plane)
 plant=ss(A,B,C,D);
 pole = eig(A);
 
 %Membuat stabil dengan menggunakan LQR
 Q=C'*C;
 R=0.01*eye(2);
 p=care(A,B,Q,R);
 K=inv(R)*B'*p;
 
 Anew=A-B*K;
 Bnew=B;
 Cnew=C;
 Dnew=D;
 polenew=eig(Anew);
 
 %==============================================
 %Membuat frequency response
tfunction=tf(plant);
bode(tfunction);
 %================================================
 %Diperoleh sistem baru yang lebih stabil
 %Perhitungan H2 norm, Hinfinity norm, L2 dan L infinity
 
 %computing H2
 plant=ss(Anew,Bnew,Cnew,Dnew);
 H2=norm(plant,2);
 H2=H2^2;
 
 Pgram=gram(Anew,Bnew);
 Qgram=gram(Anew',Cnew');
 
 %hitung h infinity
 G=pck(Anew,Bnew,Cnew,Dnew);
 hinfnorm(G,0.0001)
 linfnorm(G,0.0001)
 w=logspace(-25,25,200);
 Gf=frsp(G,w);
 [u,s,v]=vsvd(Gf);
 vplot('liv,lm',s), grid
 
 
 sysgabungan=[A B;C D];
 
 %computing H inf norm
 ninf = hinfnorm(G) %Checking H? System
 
 %computing L inf norm
n = norm(G,Inf) %Checking L? System
 
%computing L2 norm
l2 = sqrt(sum(abs(sysgabungan).^2)) %Checking L2 System
l2 = max(svd(sysgabungan)) %Checking L2 System
matL2 = norm(sysgabungan,2) %Checking L2 System