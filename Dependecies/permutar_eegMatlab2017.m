%%% Funci�n para permutar los Averages de cada sujeto para dos condiciones.

%%% Emi Fló Mayo 2017 CIBPsi

%%% Permuta aleatoriamente para cada sujeto los Averages entre si. 
%%% Es decir que para cada sujeto en cada permutacion el avg de la condicion 1 puede 
%%% intercambiarse por el avg de la condicion 2 (se shuflean las etiquetas)

%%% OJO! No admite shuflear más de 2 condiciones, tengo que arreglarlo :)

% ------------------------------------------- ------------------------
%Requiere como argumentos:

% la estructura de allsubj de fieldtrip para cada una de las condiciones:
% allsubj_trig1 y allsubj_trig2. EL ORDEN DE LOS ARGUMENTOS IMPORTA! -
% c la cantidad de sujetos
% g la cantidad de permutaciones a realizar
% cond la cantidad de condiciones

%La salida es un array de 1xg. Dentro de cada cell hay un cel de cond x
%cantidad de sujetos

%Enero 2018: Modifico s = RandStream( 'mt19937ar','seed',sum(100*clock)); por
% s = RandStream( 'mt19937ar','seed',sum(100*clock+i+j));
% otra posibilidad sería s = RandStream( 'mt19937ar','seed','shuffle');
% pero demora mucho (10 minutos para generar las permutaciones solamente).

function [ permutado_ ] = permutar_eegMatlab2017( allsubj_trig1,allsubj_trig2,c,g,cond)

tic
por_sujeto = cell(cond,1);

para_permutar = cell(cond,c); 

label_swap = cell(cond,c);

permutado_ = cell(1,g);

%la salida del script preprocesa_composition allsubj_trig1 = condici�n dos palabras 
%allsubj_trig2 = condici�n una palabra, contradictorio i know ;)

OW = allsubj_trig2; 
TW = allsubj_trig1;

for i=1:length(OW)
    para_permutar{1,i}=OW{1,i};
    para_permutar{2,i}=TW{1,i}; 
end   

for j=1:g
    
    s = RandStream( 'mt19937ar','seed','shuffle');
 
    for i=1:length(OW)

        por_sujeto{1,1}=para_permutar{1,i}; %en por_sujeto en cada iteraci�n est�n los promedios de las dos condiciones para un sujeto
        por_sujeto{2,1}=para_permutar{2,i};
        
        Perm = randperm(s,cond);
        sujeto_perm=por_sujeto(Perm);

        label_swap{1,i} = sujeto_perm{1,1}; %label_swap contiene las permutaciones para los 23 sujetos
        label_swap{2,i} = sujeto_perm{2,1}; %cada permutacion (de un total g)se almacenan en permutado_

        permutado_{1,j}=label_swap;

    end
    
    j=j+1;
end
toc


%permutado__adj=permutado_;
%save ('permutaciones_comp','permutado_')
%save ('permutaciones_adj','permutado__adj')



