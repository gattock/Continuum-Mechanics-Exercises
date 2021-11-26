clc; clear all ; lol=10^(-15); %wanna display 4 digits
a=-1.0+lol; b=3.0+lol; c=-1.0+lol; dim=3; %lol=infinitesimal error
format short %cut decimals to 4
format compact %cmd window compact out
dash='------------------------------------------------';
sigma=[a b c ; b a*c a*b ; c a*b a*b*c]
[oldevect,lamdamat]=eig(sigma)
                                                        disp(dash)
%tested that oldevect is already normalized
selfsigma=(oldevect.')*sigma*oldevect % P.' A P not P A P.'!!
                                    disp(dash)
for i=1:dim %change eig.val. into a vector
    oldeval(i)=lamdamat(i,i);
end
oldeval
oldord=['A' 'B' 'C']
                                    disp(dash)
%call sorting subroutine
[neword,eigval,eigvec]= sorting_vector_subrout(dim,oldeval,oldord,oldevect);
neword
eigval
eigvec
                                     disp(dash)
for i=1:dim %weird way to plot the sorting e.g. "C<-->N(1)"
String=[ neword(i), '<-->N(', num2str(i) ,')' ];
disp(String)
N1(i)=eigvec(i,1);
N2(i)=eigvec(i,2);
N3(i)=eigvec(i,3);
end
                                       disp(dash)
Sigma1=eigval(1)
N1
Sigma2=eigval(2)
N2
Sigma3=eigval(3)
N3
                                    disp(dash)
%check N1XN2=N3
C=cross(N1,N2)
check=C-N3
disp ('check=N1xN2-N3 should be (0 0 0)')
if check(1)+check(2)+check(3)<10^(-4)
    disp('N1 N2 N3 base is right-handed')
else
    disp('N1 N2 N3 base is left handed')
                                        disp(dash)
    disp('hence we introduce a new vector N3=-N3 and we change its eig.val')
    N3=-N3
    %eigval(3)=(-1)*eigval(3) %error, not to do! compression stays that
    check=C-N3
    disp('N1 N2 N3 new base is right-handed')
end
                                    disp(dash)
selfbase(:,1)=N1; selfbase(:,2)=N2; selfbase(:,3)=N3;
selfbase
                                    disp(dash)
%EX1B
for i=1:dim
quiver3( 0,0,0,selfbase(1,i), selfbase(2,i), selfbase(3,i) )
hold on
end
quiver3(0,0,0, 1,0,0)
quiver3(0,0,0, 0,1,0)
quiver3(0,0,0, 0,0,1)
hold off

%EX1C: I is a vector containing idexed invariants i.e. I(i)=Ii
I(1)=eigval(1)+eigval(2)+eigval(3);
I(2)=eigval(1)*eigval(2)+eigval(2)*eigval(3)+eigval(3)*eigval(1);
I(3)=eigval(1)*eigval(2)*eigval(3);
I
                                    disp(dash)
%EX1D
Smax= ( eigval(1)-eigval(3) )/2/4
 %subroutine cube_plot and stress vectors
lato=1 ;  center=[0.5 0.5 0.5];
[ciao]= cube_plot(lato,center);
norm=zeros(dim,1);
for i=1:dim
    sumsquared=0;
        for j=1:dim
        sumsquared=sumsquared+sigma(j,i)^2;
        end
    norm(i)=(sumsquared)^0.5;
end
%plot of full face stresses
quiver3(1,0.5,0.5,sigma(1,1)/norm(1,1),sigma(1,2)/norm(1,1),sigma(1,3)/norm(1,1) )
hold on
quiver3(0.5,1,0.5,sigma(2,1)/norm(2,1),sigma(2,2)/norm(2,1),sigma(2,3)/norm(2,1) )
hold on
quiver3(0.5,0.5,1, sigma(3,1)/norm(3,1),sigma(3,2)/norm(3,1),sigma(3,3)/norm(3,1) )
hold on
%plot of axial components
quiver3(1,0.5,0.5,sigma(1,1)/norm(1,1),0,0 )
hold on
quiver3(0.5,1,0.5,0,sigma(2,2)/norm(2,1),0 )
hold on
quiver3(0.5,0.5,1, 0,0,sigma(3,3)/norm(3,1) )
hold on
%plot of shearmax vector vS
vS=zeros(dim,1);
vS=[eigvec(1,1)+eigvec(1,3) eigvec(2,1)+eigvec(2,3) eigvec(3,1)+eigvec(3,3)]*Smax
quiver3(0.5, 0.5, 0.5,vS(1),vS(2),vS(3))
hold off
                                    disp(dash)
                                    disp(dash)
                                    disp('EX2;  [E]=[sigma]=GPa')
%EX2
E=110 % [GPa]
coeff=0.1
sigmaex=sigma*coeff %[GPa]
                                    disp(dash)
ni1=0.33 %metals
epsilon1=( (1+ni1)*sigmaex-ni1*eye(dim)*I(1)*coeff)/E    %ask for I(1)=tr(sigma?)
                                    disp(dash)
ni2=0.5 %incompressibles
epsilon2=( (1+ni2)*sigmaex-ni2*eye(dim)*I(1)*coeff)/E
                                    disp(dash)
                                    disp(dash)
                                    disp('EX3')
%EX3                                    
sigma
sigma_I=[1.1000 -0.2101 -3.3899 ; 
                -0.2101 0.6500 -0.5500 ; 
                -3.3899 -0.5500 1.2500]
sigma_II=[2.0000 -1.5858 -4.4142 ; 
                 -1.5858 0.5000 0.5000 ; 
                 -4.4142 -0.5000 0.5000]
                                                 disp(dash)
disp('which one is sigma rotated among I and II?')
                                    disp(dash)
[eigvectI,eigvalI]=eig(sigma_I);
[eigvectII,eigvalII]=eig(sigma_II);
eigvectI
eigvalI
                                    disp(dash)
eigvectII
eigvalII
                                    disp(dash)
 risI=0; risII=0;
for i=1:dim
    risI=risI+eigvalI(i,i)-lamdamat(i,i);
    risII=risII+eigvalII(i,i)-lamdamat(i,i);
end
risI
risII
if risI<10^(-7)
        disp('It exist an ortho-normal transformation')
        disp('that brings sigma into sigma_I')
end
if risII<10^(-7)
        disp('It exist an ortho-normal transformation')
        disp('that brings sigma into sigma_II')
end 
                                    disp(dash)     
                                    disp('HELLO WORLD!')
    