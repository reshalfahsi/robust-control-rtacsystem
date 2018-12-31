clear all;

%nilai paramter 

clear all;

M = 1.3608;

m = 0.096;

R = 0.0592;

I = 2.175*10^-4;

k = 186.3;
 


%pembagi di persamaan

h = (M+m)*(I+0.5*m*R^2);
 


%nilai di matriks A

a11 = 0;
a12 = 1;
a13 = 0;
a14 = 0;
a21 = (-k*(I+m*R^2))/h;
a22 = 0;
a23 = -(((m^2)*(R^3)*(0.555*((M+m)*I+(M*m*R^2))+m*(0.616*I+0.555*m*R^2)))/(h^2));

a24 = (1.57*I*m*R+1.11*(m^2)*(R^3))/h;

a31 = 0;
a32 = 0;
a33 = 0;
a34 = 1;
a41 = (0.707*m*R*k)/h;

a42 = 0;
a43 = 0;
a44 = -((0.785*(m^2)*(R^2))/h);
 
%nilai di matriks B
b11 = 0;
b12 = 0;
b21 = (I+m*R^2)/h;
b22 = (-0.707*m*R/h);
b31 = 0;
b32 = 0;
b41 = (-0.707*m*R/h);
b42 = (M+m)/h;
 
%matriks A,B,C, dan D
A = [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44];
B = [0 0; b21 b22; 0 0; b41 b42];
C = eye(4);
D = zeros(4,2);
%================================================================

%observable
Ob=obsv(A,C);
[n,m]=size(Ob);
unob = m-rank(Ob);
if(unob==0)
   disp('Given System is Observable.');
else
   disp('Given System is Unobservable');
end
 
%controllable
Co = ctrb(A,B);
[n,m]=size(Co);
unco=n-rank(Co);
if(unco==0)
   disp('Given System is Controllable.');
else
   disp('Given System is Uncontrollable');
end
%================================================
 

%Mengecek kestabilan (melihat apakah pole berada di half-left plane)
 
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
plant=ss(Anew,Bnew,Cnew,Dnew);



%==============================================
%Membuat frequency response
tfunction=tf(plant);
bode(tfunction);
%================================================
%Diperoleh sistem baru yang lebih stabil
%Perhitungan H2 norm, Hinfinity norm, L2 dan L infinity
 
%computing H2
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
