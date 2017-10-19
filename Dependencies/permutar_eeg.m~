%%% Funci�n para permutar los Averages de cada sujeto para dos condiciones.

%%% Emi Fló Mayo 2017 CIBPsi

%%% Permuta aleatoriamente para cada sujeto los Averages entre si. 
%%% Es decir que para cada sujeto en cada permutacion el avg de la condicion 1 puede 
%%% intercambiarse por el avg de la condicion 2 (se shuflean las etiquetas)

%%% OJO! No admite shuflear más de 2 condiciones, tengo que arreglarlo :)

% -------------------------------------------------------------------
%Requiere como argumentos:

% la estructura de allsubj de fieldtrip para cada una de las condiciones: allsubj_trig1 y allsubj_trig2
% c la cantidad de sujetos
% g la cantidad de permutaciones a realizar
% cond la cantidad de condiciones

%La salida es un array de 1xg. Dentro de cada cell hay un cel de cond x
%cantidad de sujetos

function [ permutado_ ] = permutar_eeg( allsubj_trig1,allsubj_trig2,c,g,cond)


permutado__ = cell(1,g);
por_sujeto = {};


s = RandStream.create( 'mt19937ar','seed',sum(100*clock));
RandStream.setDefaultStream(s)

para_permutar = cell(cond,c); 

%la salida del script preprocesa_composition allsubj_trig1 = condici�n dos palabras 
%allsubj_trig2 = condici�n una palabra, contradictorio i know ;)
OW = allsubj_trig2; 
TW = allsubj_trig1;

for i=1:length(OW)
    para_permutar{1,i}=OW{1,i};
    para_permutar{2,i}=TW{1,i}; 
end   

for j=1:g
    for i=1:length(OW)

        por_sujeto{1,1}=para_permutar{1,i};
        por_sujeto{2,1}=para_permutar{2,i};
        Perm = randperm(length(por_sujeto));
        sujeto_perm=por_sujeto(Perm);

        label_swap{1,i} = sujeto_perm{1,1};
        label_swap{2,i} = sujeto_perm{2,1};

        permutado_{1,j}=label_swap;
    end
    j=j+1;
end



%permutado__adj=permutado_;
%save ('permutaciones_comp','permutado_')
%save ('permutaciones_adj','permutado__adj')



