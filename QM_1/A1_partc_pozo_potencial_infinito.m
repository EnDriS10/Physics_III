clear; clc;

m=9.109E-31;
h=6.626E-34;
h_ = h / (2 * pi);
L=2;


n= (1:10);

E1= ((h_^2*pi^2)/(2*m*L));
E= E1.*n.^2;


x=(0:0.0005:L);

for i=1:length(n)
    Yn(i,:) = sqrt(2/L).*sin( (( i.*pi.*x ) ./ L) );
    PDF(i,:) = Yn(i,:).^2;
  
 figure(1)
 hold on
 subplot(6,2,i) ;
 plot(x,Yn(i,:),'b') ; % Wave f
 grid on ; % mostrar retícula

 figure(2)
 hold on
 subplot(6,2,i) ;
 plot(x,PDF(i,:),'r') ; % Plot the probability density
 grid on ; % mostrar retícula

end
   


%% Probabilidad de encontrar la particula

n_=2;
a=0.4*L;
b=0.6*L;
x_dat= x(x >= a & x <= b);

PDF_ = PDF(n_,:);
PDF_dat=PDF_(x >= a & x <= b);
[~ , Prob ] = Int(x_dat, PDF_dat, 0);

fprintf('La Probabilidad de encontar al electron entre %.2f y %.2f es de %.3f %% \n',a,b, Prob*100) ;

%% Valor esperado de la posicion

for i=1:length(n)

[~ , Valor_x(i) ] = Int(x, x.*PDF(i,:), 0);
fprintf('Para (n=%.0f) <x> = %.3fL \n',i,Valor_x(i)/L) ;

end

%%
%INTEGRAL
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end
