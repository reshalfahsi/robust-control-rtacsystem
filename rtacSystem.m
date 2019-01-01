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
%Mengecek kestabilan (melihat apakah pole berada di half-right plane)
plantlama=ss(A,B,C,D);
openloop=tf(plantlama)
pole = eig(A);
 
%Gs Gu decomposition sblm LQR
Glama=pck(A,B,C,D);
[GsLAMA,GuLAMA]=sdecomp(Glama);
 
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
 
%================================================
%Diperoleh sistem baru yang lebih stabil
%Perhitungan H2 norm, Hinfinity norm, L2 dan L infinity
 
G=pck(Anew,Bnew,Cnew,Dnew);
plantbaru=ss(Anew,Bnew,Cnew,Dnew);
closedloop=tf(plantbaru);

%Gs Gu decomposition stlh LQR
[Gs,Gu]=sdecomp(G)

%Observability Gramian Matrix
Qgram=gram(cjt(Anew),cjt(Cnew));

%computing H2 and L2
H2=norm(plantbaru,2)
% atau
H2=h2norm(G) 
L2=sqrt(trace(cjt(Bnew)*Qgram*Bnew))
 
%hitung h infinity dan l infinity
Hinf = hinfnorm(G,0.0001)
Linf = linfnorm(G,0.0001)
 
%% Internal Stability
syms s real% Symbolic variable 's'

% Open Loop System
P = minreal(openloop) %Reduce to minimum order

for i=1:1:4 % Convert transfer function to symbolic
 for j=1:1:2
     symP(i,j) = tf2sym(P,i,j);
     symG_cl(i,j) = tf2sym(closedloop,i,j);
 end
end

%% Dari persamaan mencari controller K jika diketahui system closed loop G_cl dan plant P
%% K = [(G(s)^-1)P(s) - I] P(s)^-1
symK = simplify((pinv(symG_cl)*symP*pinv(symP)- pinv(symP)));
symK_hat = -symK;
for i=1:1:2 % Convert symbolic back to transfer function
 for j=1:1:4
 K_hat(i,j) = syms2tf(symK_hat(i,j));
 end
end
K_hat

%% RHP poles test
symPK_hat = symP*symK_hat; % P(s)K_hat(s)
for i=1:1:4 % Convert symbolic back to transfer function
 for j=1:1:4
  PK_hat(i,j) = syms2tf(symPK_hat(i,j));
 end
end

%% (I-P*K_hat)^-1
symIStfmatrix = simplify(inv(eye(size(symP*symK_hat)))-symP*symK_hat)
for i=1:1:4 % Convert symbolic back to transfer function
 for j=1:1:4
 IStfmatrix(i,j) = syms2tf(symIStfmatrix(i,j));
 end
end

if(isstable(IStfmatrix)==0)
 disp('Not Internally Stable')
else
 disp('Internally Stable')
end
 
%% Coprime Factorization
syms s real
F = -lqr(A, B, Q, R); %gain F sehingga A+BF stabil
%% M(s)
AM = A+B*F;
BM = B;
CM = F'*F;
DM = eye(size(D));
ssM = ss(AM,BM,CM,DM); % state space of system M
tfM = tf(ssM) % transfer function of system M
 
%% N(s)
AN = A+B*F;
BN = B;
CN = C+D*F;
DN = D;
ssN = ss(AN,BN,CN,DN); % state space of system N
tfN = tf(ssN) % transfer function of system N
 
%Convert to Symbolic
for i=1:1:4
 for j=1:1:2
     symN(i,j) = tf2sym(tfN,i,j);
     symM(i,j) = tf2sym(tfM,i,j);
 end
end

% diuji apakah P=N*M^-1
symPcop = symN*pinv(symM)

%Convert symbolic back to transfer function
for i=1:1:4
 for j=1:1:4
 Pcop(i,j) = syms2tf(symPcop(i,j));
 end
end
Pcop