function [orden,rho,FPE,AIC,AICm,a,sigma,ref]=arord(y,n,burg,criterio)

%ARORD Hace una estimación del orden del modelo y calcula los
%  parametros del modelo AR para la secuencia 'y' utilizando el
%  metodo recursivo de Burg.
%
%	   [ORDEN,RHO,FPE,AIC,AICm,A,SIGMA,REF]=ARORD(Y,N,BURG,CRITERIO)
%
%  ORDEN: Orden seleccionado.
%
%  RHO: Potencias del error de predicción estimadas para cada orden
%      (RHO(N+1)=SIGMA).
%
%  FPE: Error de predicción final (útil para la selección del orden).
%
%  AIC: Criterio de Información de Akaike (útil para la selección
%       del orden).
%
%  AICm: Criterio de Información de Akaike modificado (útil para la
%         selección del orden).
%
%  A: Vector de parametros del modelo AR de orden N.
%
%  SIGMA: Potencia del error de predicción estimada (varianza del error)
%
%  REF: Coeficientes de reflexión estimados.
%
%  Y: Secuencia de salida del modelo .
%
%  N: Orden que se considera máximo para el modelo AR.
%
%
%  BURG: Permite indicar si se quiere utilizar la media armónica
%      de la suma de los cuadrados de los errores de predicción
%      (método de Burg: BURG=1), o bien la geométrica (BURG=0).
%
%  CRITERIO Criterio de selección del orden utilizado:
%       CRITERIO=0 -> Utilizando RHO.
%       CRITERIO=1 -> Utilizando FPE.
%       CRITERIO=2 -> Utilizando AIC.
%       CRITERIO=3 -> Utilizando AICm.


%---------------------------------------------------------------
% Eliminamos el valor medio de la secuencia:
y=detrend(y,0);
%---------------------------------------------------------------
% Por defecto se utiliza el metodo de Burg (media armonica) y el
% criterio de Akaike (AIC)
if nargin<3,
  burg=1;
end
if nargin<4,
  criterio=2;
end

%---------------------------------------------------------------
% Orientamos 'y' como un vector columna
y=y(:);
L=length(y);

%---------------------------------------------------------------
% Calculamos los coeficientes AR utilizando el algoritmo de Burg:

ef=y;
eb=y;
rho(1)=y'*y/L;

for p=1:n
   nef=ef(p+1:L)'*ef(p+1:L);
   neb=eb(p:L-1)'*eb(p:L-1);
   if burg,
     den=(nef+neb)/2;
   else
     den=sqrt(nef*neb);
   end
   ref(p)=(-eb(p:L-1)'*ef(p+1:L))/den;
   efold=ef;
   ef(2:L)=ef(2:L)+ref(p)*eb(1:L-1);
   eb(2:L)=eb(1:L-1)+conj(ref(p))*efold(2:L);
   a(p)=ref(p);
   a(1:p-1)=a(1:p-1)+ref(p)*conj(a(p-1:-1:1));
   rho(p+1)=rho(p)*(1-ref(p)*ref(p));
end

a=[1 a];

sigma=rho(n+1);

FPE=zeros(1,n+1);
AIC=zeros(1,n+1);
AICm=zeros(1,n+1);

for i=1:n+1,
  FPE(i)=rho(i)*(1+i/L)/(1-i/L);
  AIC(i)=L*log(rho(i))+2*i;
  AICm(i)=L*log(rho(i))+i*log(L);
end

if criterio==0,
  rholim=0.005*(max(rho)-min(rho))+min(rho);
  [valormin,posmin]=min(abs(rho-rholim));
  orden=posmin-1; 
elseif criterio==1,
    [valormin,posmin]=min(FPE);
    orden=posmin-1;    
elseif criterio==2,
    [valormin,posmin]=min(AIC);
    orden=posmin-1;
elseif criterio==3,
    [valormin,posmin]=min(AICm);
    orden=posmin-1;
end

