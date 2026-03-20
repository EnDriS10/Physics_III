clear; clc;

%% DATOS
m=9.109E-31;
h=6.626E-34;
h_ = h / (2 * pi);

V0=1;

dx = 1E-4;

var = 4; % var = (m*V0*L^2)/(2*h_^2)^0.5

%% CLASE 1

p1= @(x) x.*tan(x) - ( (var^2) - x.^2 ).^0.5;
x=0:dx:2*pi;
[ep1]=mi_zeros_autom(x,p1(x));

ep1=ep1( abs(cos(ep1)) > dx  & ((var^2) - ep1.^2 ) > 0 );
E1 = V0.*(ep1./var).^2;

alf1= sqrt( (2.*m.*(V0-E1))./(h_^2));
k1= sqrt( (2.*m.*E1)./(h_^2));

%% CLASE 2
 
p2= @(x) x.* cot(x) + ( (var^2) - x.^2 ).^0.5;
[ep2]=mi_zeros_autom(x,p2(x));

ep2=ep2( abs(sin(ep2)) > dx  & ((var^2) - ep2.^2 ) > 0 );
E2 = V0.*(ep2./var).^2;

alf2= sqrt( (2.*m.*(V0-E2))./(h_^2));
k2= sqrt( (2.*m.*E2)./(h_^2));

%% Energias y Epsilons

ep =sort([ep1, ep2]);
E =sort([E1, E2]);


%% Ec Onda Sin Normalizar

L= sqrt( (2*h_^2)*(var^2) / (m*V0) ); % L en funcion de V0

x_1= -10*L:L*dx:-L/2;
x_2= -L/2:L*dx:L/2;
x_3= L/2:L*dx:10*L;

% Onda 1
y1_1= cos(k1(1)*L/2)*exp(alf1(1)*L/2).*exp(alf1(1).*x_1);
y1_2= cos(k1(1).*x_2);
y1_3= cos(k1(1)*L/2)*exp(alf1(1)*L/2).*exp(-alf1(1).*x_3);

%PDF 1
pdf1_1 = y1_1.^2;
pdf1_2 = y1_2.^2;
pdf1_3 = y1_3.^2;

% Onda 2
y2_1= -sin(k2(1)*L/2)*exp(alf2(1)*L/2).*exp(alf2(1).*x_1);
y2_2= sin(k2(1).*x_2);
y2_3= sin(k2(1)*L/2)*exp(alf2(1)*L/2).*exp(-alf2(1).*x_3);

%PDF 2
pdf2_1 = y2_1.^2;
pdf2_2 = y2_2.^2;
pdf2_3 = y2_3.^2;


% Onda 3
y3_1= cos(k1(2)*L/2)*exp(alf1(2)*L/2).*exp(alf1(2).*x_1);
y3_2= cos(k1(2).*x_2);
y3_3= cos(k1(2)*L/2)*exp(alf1(2)*L/2).*exp(-alf1(2).*x_3);

%PDF 3
pdf3_1 = y3_1.^2;
pdf3_2 = y3_2.^2;
pdf3_3 = y3_3.^2;

%% Ec Onda Normalizada

[~, C1_1 ] = Int(x_1,pdf1_1, 0);
[~, C1_2 ] = Int(x_2,pdf1_2, 0);
[~, C1_3 ] = Int(x_3,pdf1_3, 0);
B1 =1/sqrt( C1_1 + C1_2 + C1_3);

[~, C2_1 ] = Int(x_1,pdf2_1, 0);
[~, C2_2 ] = Int(x_2,pdf2_2, 0);
[~, C2_3 ] = Int(x_3,pdf2_3, 0);
A1 = 1/sqrt( C2_1 +C2_2 +C2_3 );

[~, C3_1 ] = Int(x_1,pdf3_1, 0);
[~, C3_2 ] = Int(x_2,pdf3_2, 0);
[~, C3_3 ] = Int(x_3,pdf3_3, 0);
B2 = 1/sqrt(C3_1 +C3_2 +C3_3 );

% Normalizing the wave functions
y1_1 = B1 * y1_1;
y1_2 = B1 * y1_2;
y1_3 = B1 * y1_3;

y2_1 = A1 * y2_1;
y2_2 = A1 * y2_2;
y2_3 = A1 * y2_3;

y3_1 = B2 * y3_1;
y3_2 = B2 * y3_2;
y3_3 = B2 * y3_3;
 
%% Plot

figure(1)
hold on
h1=plot(x_1./L, y1_1.*sqrt(L), 'r',x_2./L, y1_2.*sqrt(L),'r',x_3./L, y1_3.*sqrt(L),'r' );
h2=plot(x_1./L, y2_1.*sqrt(L), 'b',x_2./L, y2_2.*sqrt(L),'b',x_3./L, y2_3.*sqrt(L),'b' );
h3=plot(x_1./L, y3_1.*sqrt(L), 'g',x_2./L, y3_2.*sqrt(L),'g',x_3./L, y3_3.*sqrt(L),'g' );
xl1=xline(-0.5, 'k--', 'LineWidth', 1.5);
xl2=xline( 0.5, 'k--', 'LineWidth', 1.5);
grid on

xlabel('Posición (x/L)');
ylabel('Wave Function \psi \cdot \surdL');
title('Wave Functions');
legend( ...
    [h1(1), h2(1), h3(1), xl1], ...
    {['\psi_1 (par, E_1/V_0 = '  + string(round(E1(1)/V0, 3)) + ')'], ...
     ['\psi_2 (impar, E_2/V_0 = ' + string(round(E2(1)/V0, 3)) + ')'], ...
     ['\psi_3 (par, E_3/V_0 = '  + string(round(E1(2)/V0, 3)) + ')'], ...
     'Paredes de la caja'}, ...
    'Location', 'northeast');

figure(2)

% --- Onda 1 ---
subplot(3,2,1)
hold on
plot(x_1./L, y1_1.*sqrt(L), 'r')
plot(x_2./L, y1_2.*sqrt(L), 'r')
plot(x_3./L, y1_3.*sqrt(L), 'r')
xline(-0.5, 'k--'); xline(0.5, 'k--');
title('\psi_1 (par)')
xlabel('x/L'); ylabel('\psi \cdot \surdL')
grid on; hold off

subplot(3,2,2)
hold on
plot(x_1./L, (y1_1.*sqrt(L)).^2, 'r')
plot(x_2./L, (y1_2.*sqrt(L)).^2, 'r')
plot(x_3./L, (y1_3.*sqrt(L)).^2, 'r')
xline(-0.5, 'k--'); xline(0.5, 'k--');
title('|\psi_1|^2 (par)')
xlabel('x/L'); ylabel('|\psi|^2 \cdot L')
grid on; hold off

% --- Onda 2 ---
subplot(3,2,3)
hold on
plot(x_1./L, y2_1.*sqrt(L), 'b')
plot(x_2./L, y2_2.*sqrt(L), 'b')
plot(x_3./L, y2_3.*sqrt(L), 'b')
xline(-0.5, 'k--'); xline(0.5, 'k--');
title('\psi_2 (impar)')
xlabel('x/L'); ylabel('\psi \cdot \surdL')
grid on; hold off

subplot(3,2,4)
hold on
plot(x_1./L, (y2_1.*sqrt(L)).^2, 'b')
plot(x_2./L, (y2_2.*sqrt(L)).^2, 'b')
plot(x_3./L, (y2_3.*sqrt(L)).^2, 'b')
xline(-0.5, 'k--'); xline(0.5, 'k--');
title('|\psi_2|^2 (impar)')
xlabel('x/L'); ylabel('|\psi|^2 \cdot L')
grid on; hold off

% --- Onda 3 ---
subplot(3,2,5)
hold on
plot(x_1./L, y3_1.*sqrt(L), 'g')
plot(x_2./L, y3_2.*sqrt(L), 'g')
plot(x_3./L, y3_3.*sqrt(L), 'g')
xline(-0.5, 'k--'); xline(0.5, 'k--');
title('\psi_3 (par)')
xlabel('x/L'); ylabel('\psi \cdot \surdL')
grid on; hold off

subplot(3,2,6)
hold on
plot(x_1./L, (y3_1.*sqrt(L)).^2, 'g')
plot(x_2./L, (y3_2.*sqrt(L)).^2, 'g')
plot(x_3./L, (y3_3.*sqrt(L)).^2, 'g')
xline(-0.5, 'k--'); xline(0.5, 'k--');
title('|\psi_3|^2 (par)')
xlabel('x/L'); ylabel('|\psi|^2 \cdot L')
grid on; hold off

sgtitle('Funciones de onda y densidades de probabilidad - Pozo finito (var = 4)')


figure(3)
hold on

% --- Potencial en forma de caja (normalizado por V0) ---
x_pot = [-3, -0.5, -0.5,  0.5,  0.5,  3];
V_pot = [ 1,   1,    0,    0,    1,    1];
plot(x_pot, V_pot, 'k', 'LineWidth', 2)

% --- Líneas de energía dentro del pozo ---
colores = {'r', 'b', 'g', 'm'};
for i = 1:length(E)
    Ei = E(i)/V0;
    col = colores{mod(i-1, length(colores)) + 1};
    plot([-0.5, 0.5], [Ei, Ei], '--', 'Color', col, 'LineWidth', 1.5)
    text(0.52, Ei, ['E_' num2str(i) '/V_0 = ' num2str(round(Ei,3))], ...
        'Color', col, 'FontSize', 9, 'VerticalAlignment', 'middle')
end

% --- Formato ---
xlim([-3, 3])
ylim([-0.05, 1.2])
xlabel('x/L')
ylabel('V / V_0')
title('Pozo finito: potencial y espectro de energía  ')
grid on
hold off
%%  FUNCIONES

 function [zeros]=mi_zeros_autom(x,y)

ff=y(1:end-1).*y(2:end) ; % atención, tiene un elemento menos
 x_mid=(x(1:end-1)+x(2:end))/2 ; % puntos medios de los intervalos
 b_negativos=ff<0 ;      % generar vector lógico
 zeros=x_mid(b_negativos); % extraer b_negativos de x_mid

 end

%INTEGRAL
function [Int_Indef, Int_def ] = Int(x,y, y0)

Int_Indef =zeros(size(y)) ;
areas=diff(x).*(y(1:end-1)+y(2:end))/2 ;
Int_def=sum(areas) ;
Int_Indef=[0, cumsum(areas)];

Int_def=y0 + Int_def ;
Int_Indef=y0 + Int_Indef;

end
