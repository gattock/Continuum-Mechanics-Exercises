clc; clear all; %ver3
format shortG %cut decimals as needed
format compact %cmd window compact out
%axis equal
a=-1; b=1; c=0; d=3; e=3; f=0; g=0; h=0; l=1; 
aa=140; bb=160; kk=300; %alpha,beta,kappa
dash="--------------------------------------------------------------------------";
T=[aa bb kk ; bb aa+bb bb+kk ; kk bb+kk aa+bb+kk]
disp("Ex1: compute deformation gradient tensor")
F=[ a b c; d e f ; g h l]
disp(dash)
disp("Ex2:")%ex2:compute 2nd piola-kirchhoff stress tensor S
invF=inv(F)
J=det(F)
disp("------------------------")
disp("First PK stress tensor")
T0=J * T * invF' %first piola-kirch stress tens
disp("Second PK stress tensor")
S=invF *T0 %second p-k 

disp(dash) %ex3 Von Mises
disp('Ex3:Von Mises stress computation')
[eigvec,eigval]=eig(T)
trT=trace(T)
T_hyd=zeros(3,3);
for i=1:3
T_hyd(i,i)=trT/3;
end
T_hyd
T_dev=T-T_hyd

temp=0;
for i=1:3 %computation of sigma prime
    for j=1:3
        temp=temp+T_dev(i,j)^2;
    end
end
sigma_vm=(3/2*temp)^0.5

disp(dash) %ex4 
disp("Ex4: Strain tensor Epsilon computation by Hooke's law")
E=70 %[GPa]
nu=0.3
Eps=[(1+nu)*T-nu*eye(3)*trT]/E
disp("Show the relation between deviators Tˆ_ij = 2 G Epsˆ_ij")
T_dev
for i=1:3
Eps_hyd(i,i)=trace(Eps)/3;
end
Eps_hyd
Eps_dev=Eps-Eps_hyd
G=E/2/(1+nu)
disp("2 G Eps_dev=")
2*G*Eps_dev
TestVar=det(T_dev-2*G*Eps_dev)
if TestVar<10^(-9)
   disp("The relation between deviators is verified")
end

disp(dash) %ex5
disp("Ex5: Stiffness tensor computation")

%ijkq=2323
I=eye(3);
for i=1:3
    for j=1:3
        for k=1:3
            for q=1:3
C(i,j,k,q)=E/(1+nu)*(0.5*( I(i,k)*I(j,q)+I(j,k)*I(i,q))+(nu/(1-2*nu)*I(i,j)*I(k,q)));
            end
        end
    end
end
C
disp(dash)
disp("C(2,3,2,3) is equal to")
C(2,3,2,3)
disp(dash)
disp("Hello World")



