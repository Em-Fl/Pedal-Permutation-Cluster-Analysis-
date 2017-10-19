 
% Emi Fló junio 2017 @CIBPsi
%% -- Script para Permutation Cluster Analisis siguiendo el analisis realizado por (Bemis&Pylkaneen 2011). 

%%% Pylka usa CPA :
%       -en ventana de tiempo 0-500 ms
%       -10.000 permutaciones
%       -Identifican cluster de tiempo adyacentes en los que haya una
%       interacción entre tarea y cantidad de palabras mediante un anova
%       2x2. Por lo menos 10 puntos de tiempo contiguos con  p<0.30.
%       - Para cada punto de tiempo de cada cluster ttest dentro de cada
%       tarea (factores una palabra vs dos palabras)
%       -Estadístico = SumaTvaluesComposicion -SumaTvaluesLista (Aumento
%       en la actividad durante la condicion dos palabras en comparación
%       con una palabra en la tarea de composición. Ninguna diferencia de
%       actividad entre las condiciones en la tarea lista.

%   HACE OTRA PERMUTACION ??

%%% Este script usa las funciones:

%       -permutar_eeg
%       -ttestPareadoParaValidarConFieldtrip
%       -clusterTiempo
%       -indicesClusters
%       -checkCluster2
%       -clusterTiempo_Espacio
%       -ttestClusterTiempoEspacio.

%%% Estas funciones se pueden usar para otros estadísticos
%%% siempre y cuando se parta de una matriz de pvalues (o lo que fuera) de ce x sf 
%%% (cantidad de electrodos x frecuencia de muestreo)

%%% Necesito array de vecinos. 

%%% Necesito data de sujeto como average, en fieldtrip salida allsubj_trig%

%% 0) Preparo data

% cargo data con estructura allsubj_trig de fieldtrip. Las funciones están
% diseñadas para correr sobre la estructura de las permutaciones por lo que
% la data tendrá la misma estructura. Array de 2 x cantidad de sujetos.

clc
clear

load ERP_23_subjects_adj_060517_1357
load ERP_23_subjects_comp_041617_1613
load ERP_23_subjects_sust_071417_1743

Comp2W = {};
Comp1W = {};

ListaA2W = {};
ListaA1W = {};

ListaS2W = {};
ListaS1W = {};

for i=1:length(allsubj_trig1)
    Comp2W{1,i} = allsubj_trig1{1,i}.avg([1:64],[1:512]);
    Comp1W{1,i} = allsubj_trig2{1,i}.avg([1:64],[1:512]);
    ListaA2W{1,i} = allsubj_trig_two_word_lista_adj{1,i}.avg([1:64],[1:512]);
    ListaA1W{1,i} = allsubj_trig_one_word_lista_adj{1,i}.avg([1:64],[1:512]);
    ListaS2W{1,i} = allsubj_trig_two_word_lista_sust{1,i}.avg([1:64],[1:512]);
    ListaS1W{1,i} = allsubj_trig_one_word_lista_sust{1,i}.avg([1:64],[1:512]);
end

%Guardo data real en el mismo formato que las permutaciones
DataRealComp = {};

DataRealComp(1,1:23)= Comp2W(1:23);
DataRealComp(2,1:23) = Comp1W(1:23);


DataRealListaAdjetivo = {};

DataRealListaAdjetivo(1,1:23) = ListaA2W(1:23);
DataRealListaAdjetivo(2,1:23) = ListaA1W(1:23);


DataRealListaSustantivo = {};

DataRealListaSustantivo(1,1:23) = ListaS2W(1:23);
DataRealListaSustantivo(2,1:23) = ListaS1W(1:23);
% save('DataRealComp_257puntosTiempo','DataRealComp')
% save('DataRealListaAdjetivo_257puntosTiempo','DataRealListaAdjetivo')

%% 1) Permutaciones

%%% Usa permutar_eeg: Permuta aleatoriamente para cada sujeto los Averages entre si. 
%%% Es decir que para cada sujeto en cada permutacion el avg de la condicion 1 puede 
%%% intercambiarse por el avg de la condicion 2 (se shuflean las etiquetas)

% -------------------------------------------------------------------
%Requiere como argumentos:

% la estructura de allsubj de fieldtrip para cada una de las condiciones: allsubj_trig1 y allsubj_trig2
% c la cantidad de sujetos
% g la cantidad de permutaciones a realizar
% cond la cantidad de condiciones 

%La salida es un array de 1xg. Dentro de cada cell hay un cel de cond x cantidad de sujetos. 

c=23; %cantidad de sujetos
g=2500; %cantidad de permutaciones a realizar
cond = 2;

permutado_comp = permutar_eegMatlab2017(Comp2W,Comp1W,c,g,cond);
permutado_adj = permutar_eegMatlab2017(ListaA2W,ListaA1W,c,g,cond);
% permutado_sust = permutar_eegMatlab2017(ListaS2W,ListaS1W,c,g,cond);
%save ('permutaciones_adj','permutado')
% tic
% save ('permutaciones_adj_5000','permutado_adj')
% display('pronto')
% toc


%% ANALISIS POSTA


tclusters = {}; % Array en la que guardo los estadísticos obtenidos de cada permutacion
tclusters_comp_adj = []; % Matriz en la que guardo los t de cada cluster
sf=512; %frecuencia de muestreo
ce=64; %cant de electrodos
c=23;%cant de sujetos

t0=tic;

for numeroPerm=1:1%length(permutado_comp)
    
    tic
    
    display(numeroPerm)
    
   
    %%% ---------- 2)  ANOVA 2x2 para cada punto de tiempo
    
%     para correr sobre los datos reales
%     permutacion_comp=DataRealComp;
%     permutacion_adj=DataRealListaAdjetivo;
% % % 
%     permutacion_comp = permutado_comp{numeroPerm};
%     permutacion_adj = permutado_adj{numeroPerm};
%    
    
    FACTNAMES = {'Tarea', 'CantPalabras'};

    [ pvalues_ ] = anova2x2Permutacion_beta(c,sf,ce,FACTNAMES,permutacion_comp,permutacion_adj);

    %%% ---------  3) Encuentro los clusters de tiempo-----------------------------
    
    %%% Usa funcion clusterTiempo
    
    %%% Hay que definir el valor de p y la cantidad mínima de puntos
    %%% consecutivos de tiempo con un valor de p menor al establecido para
    %%% conformar un cluster temporal
    
    threshold = 0.10; % valor p umbral de significancia
    t = 10; %cantidad de puntos de tiempo mínimos que deben conformar un cluster
    % pvalues_ matriz de ce x sf obtenida del análisis estadístico
    umbralElectrodos = 3; %mínimo de electrodos que deben formar parte del cluster para c/punto de tiempo
    
    [indices_tiempo cluster_tiempo] = clusterTiempo(ce,sf,threshold,t,pvalues_);
    
    %%% -------- 4) Defino listaElectroClusters matriz nClusters x 3 (E, a,b)----------
    %%% Estructura que preciso para funcion checkCluster
    
    % E=Electrodo donde está el cluster
    % a=indice de tiempo donde comienza cluster
    % b=indice de tiempo donde termina cluster
    
    %Requiere indices_tiempo - array 1x 64 que contiene para cada electrodo un cell array con vectores con
    %los indices de tiempo de los clusters
    
    %%% ---------- 5) Defino listaClusters -----------------------------------
    %%% cell array de electrodos, cada cell lista de indices de principio y fin de cluster (a,b).
    %%% Estructura que preciso para funcion checkCluster
    
    [ listaElectroClusters listaClusters ] = indicesClusters(indices_tiempo, ce);
 
    %%% ---------- 6) Encuentro clusters con funcióncheckCluster ----------------
    
    load vecinos_biosemi64
    
    ptable = cluster_tiempo; % Renombro matriz ce x sf que tiene 1s donde el valor de p es menor al umbral establecido
    
    [clustSoFar,clusterMatrix]=checkCluster2(ptable, vecinos,{},0,{},listaElectroClusters, listaClusters,[],0,umbralElectrodos);

    
    %%% ---------- 7) Genero una matriz ce x sf para cada cluster
  
    [ clusterTiempoEspacio ] = clusterTiempo_Espacio(clustSoFar, ce, sf);
    
    %%% ---------- 8) Ttest para cada tarea, condicion dos palabras vs una
    % palabra para cada punto de tiempo de cada cluster
    
    alpha = 0.05;
    
    [ clusterTiempoEspacioTT_comp ] = ttestClusterTiempoEspacio(permutacion_comp, clusterTiempoEspacio, c, ce, sf, alpha);
    
    [ clusterTiempoEspacioTT_adj ] = ttestClusterTiempoEspacio(permutacion_adj, clusterTiempoEspacio, c, ce, sf, alpha);

    clusters_comp_adj= cell(1,length(clusterTiempoEspacioTT_comp));
    
    for i=1:length(clusterTiempoEspacioTT_comp)
        clusters_comp_adj{1,i}= abs(clusterTiempoEspacioTT_comp{1,i})-abs(clusterTiempoEspacioTT_adj{1,i});
    end
    
%     save('clusters_compVSlistaAdj_411_3e_t010_010_005','clusters_comp_adj')
    
    %%% --------9) Construyo estadístico. Suma de tvalues de tarea 
    % composicion menos suma de tvalues de tarealista----------
    
    diferencia=zeros(1,length(clusterTiempoEspacioTT_comp));
    
    for i=1:length(clusterTiempoEspacioTT_comp)
        diferencia(1,i)=sum(abs(clusterTiempoEspacioTT_comp{1,i}(:)))-sum(abs(clusterTiempoEspacioTT_adj{1,i}(:)));
    end
    
    tclusters_comp_adj = diferencia;
%     tclusters{numeroPerm} = diferencia;
    
  toc 
  
end

t=toc(t0)

% load handel; sound(y,Fs) 
% 
% save ('tclusters_compVSlistaAdj_411_3e_t010_010_005','tclusters_comp_adj')


% save ('tclusters_2500p_compVSlistaAdj_411_3e_t010_010_005','tclusters')


%% Testeo si hay clusters significativos en mi data
clc
clear

load tclusters_2500p_compVSlistaAdj_411_3e_t010_010_005.mat
load tclusters_compVSlistaAdj_411_3e_t010_010_005.mat
load clusters_compVSlistaAdj_411_3e_t010_010_005.mat

tic
distribucion_clusters = [];

h=1;
for i=1:length(tclusters)
    for j=1:length(tclusters{i})
        distribucion_clusters(1,h)=tclusters{i}(j);
        h=h+1;
    end
end
toc

hist(distribucion_clusters)

Media = mean(distribucion_clusters);
St = std(distribucion_clusters);

cluster_significativo=[]

j=0;
for i=1:length(tclusters_comp_adj)
    [h,p,ci,zval]=ztest(tclusters_comp_adj(i),Media,St)
    if h==1
        cluster_significativo(j+1,1:2) = [i p];
        j=j+1;
    end
end

