
% Emi Fló junio 2017 @CIBPsi

%% -- Script para validar Permutation Cluster Analisis


%%% permutar_eeg, ttestPareadoParaValidarConFieldtrip, clusterTiempo, indicesClusters, checkCluster2
%%% clusterTiempo_Espacio y ttestClusterTiempoEspacio. 
%%% Estas funciones se pueden usar para otros estadísticos
%%% siempre y cuando se parta de una matriz de pvalues (o lo que fuera) de ce x sf 
%%% (cantidad de electrodos x frecuencia de muestreo)

%%% Necesito array de vecinos. 

%%% Necesito data de sujeto como average, en fieldtrip salida allsubj_trig%

%% Preparo data

% cargo data con estructura allsubj_trig de fieldtrip. Las funciones están
% diseñadas para correr sobre la estructura de las permutaciones por lo que
% la data tendrá la misma estructura. Array de 2 x cantidad de sujetos.
 clc
 clear
 
load ERP_23_subjects_adj_060517_1357
load ERP_23_subjects_comp_041617_1613

Comp2W = {};
Comp1W = {};

ListaA2W = {};
ListaA1W = {};

for i=1:length(allsubj_trig1)
    Comp2W{1,i} = allsubj_trig1{1,i}.avg([1:64],[102:512]);
    Comp1W{1,i} = allsubj_trig2{1,i}.avg([1:64],[102:512]);
    ListaA2W{1,i} = allsubj_trig_two_word_lista_adj{1,i}.avg([1:64],[102:512]);
    ListaA1W{1,i} = allsubj_trig_one_word_lista_adj{1,i}.avg([1:64],[102:512]);
end

%Guardo data real en el mismo formato que las permutaciones
DataRealComp = {};

DataRealComp(1,1:23)= Comp2W(1:23);
DataRealComp(2,1:23) = Comp1W(1:23);


DataRealListaAdjetivo = {};

DataRealListaAdjetivo(1,1:23) = ListaA2W(1:23);
DataRealListaAdjetivo(2,1:23) = ListaA1W(1:23);


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
g=5000; %cantidad de permutaciones a realizar
cond= 2;


% permutado = permutar_eegMatlab2017(Comp2W,Comp1W,c,g,cond);
permutado = permutar_eegMatlab2017(ListaA2W,ListaA1W,c,g,cond);


%save ('permutaciones_adj','permutado')
%save ('permutaciones_adj','permutado_adj')


%% 


tclusters_comp={};
tclusters={};
sf=411; %frecuencia de muestreo
ce=64; %cant de electrodos
c=23;%cant de sujetos
alpha = 0.05;

t0=tic

for numeroPerm=1:1%length(permutado)
    
    tic
    
    display(numeroPerm)
    
    %%% ---------- 2) Ttest pareado para cada punto de tiempo para composicion------------
    
    %%% Usa función ttestPareadoParaValidarConFieldtrip
    permutacion = DataRealListaAdjetivo;
    
%     permutacion = permutado{numeroPerm};

    [pvalues_] = ttestPareadoParaValidarConFieldtrip(permutacion, c, ce, sf,alpha);

    %%% ---------  3) Encuentro los clusters de tiempo-----------------------------
    
    %%% Usa funcion clusterTiempo
    
    %%% Hay que definir el valor de p y la cantidad mínima de puntos
    %%% consecutivos de tiempo con un valor de p menor al establecido para
    %%% conformar un cluster temporal
    
    threshold = 0.10; % valor p umbral de significancia
    t = 20; %cantidad de puntos de tiempo mínimos que deben conformar un cluster
    % pvalues_ matriz de ce x sf obtenida del análisis estadístico
    
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
    
    [clustSoFar,clusterMatrix]=checkCluster2(ptable, vecinos,{},0,{},listaElectroClusters, listaClusters,[]);

    %%% ---------- 7) Genero una matriz ce x sf para cada cluster
  
    [ clusterTiempoEspacio ] = clusterTiempo_Espacio(clustSoFar, ce, sf);
    
    %%% ---------- 8) Para validar con fieldtrip. Ttest para cluster de
    %%% tiempo y espacio
    
    [ clusterTiempoEspacioTT ] = ttestClusterTiempoEspacio( permutacion, clusterTiempoEspacio, c, ce, sf, alpha);
    
%     % 
    clusters_adj= clusterTiempoEspacioTT;
    
    save('clusters_adj_followUpPerm_411puntosTiempo_t20_10_005','clusters_adj')
%     
    %%% --------9) Me quedo con el maxcluster, aquel que tenga  suma de los tvalues sea mayor----------
    
    
     suma=zeros(1,length(clusterTiempoEspacioTT));
    
    for i=1:length(clusterTiempoEspacioTT)
        suma(1,i)=sum(clusterTiempoEspacioTT{1,i}(:));
    end
    
   tclusters_adj = suma; 
%    tclusters{numeroPerm} = suma;
    
  toc 
  
end


% load handel; sound(y,Fs) 

save ('tclusters_followUpPermAdjetivo_411puntosTiempo_t20_010_005','tclusters_adj')
% save ('tclusters_5000p_followUpPermAdjetivo_411puntosTiempo_t20_010_005','tclusters')

t=toc(t0)

%% Testeo si hay clusters significativos en mi data
clc
clear

load tclusters_followUpPermAdjetivo_411puntosTiempo_t20_010_005.mat
load tclusters_5000p_followUpPermAdjetivo_411puntosTiempo_t20_010_005.mat
load clusters_adj_followUpPerm_411puntosTiempo_t20_10_005.mat

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

cluster_significativo=[];

j=0;
for i=1:length(tclusters_adj)
    [h,p,ci,zval]=ztest(tclusters_adj(i),Media,St)
    if h==1
        cluster_significativo(j+1,1:2) = [i p];
        j=j+1;
    end
end


