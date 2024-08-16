function [RRf,tRRf,RRout,tRRout]=filterRR(RR,tRR,numciclos,porcentaje)

if nargin<4,porcentaje=0.2;end
if nargin<3,numciclos=1;end

elimina=iteraRR(RR,numciclos,porcentaje);

sumelimina=0;
i=1;
limitei=20;
while sum(elimina)>sumelimina & i<limitei,
    sumelimina=sum(elimina);
    elimina(~elimina)=iteraRR(RR(~elimina));
    i=i+1;
end

RRf=RR(~elimina);
tRRf=tRR(~elimina);
RRout=RR(elimina);
tRRout=tRR(elimina);


function elimina=iteraRR(RR,numciclos,porcentaje)

if nargin<3,porcentaje=0.1;end
if nargin<2,numciclos=5;end

RRm=conv(RR,[ones(1,numciclos) 0 ones(1,numciclos)]/2/numciclos);
RRm=RRm(numciclos+1:end-numciclos);
RRm([1:numciclos end-numciclos+1:end])=RR([1:numciclos end-numciclos+1:end]);
elimina=abs(RR-RRm)>porcentaje*RRm;