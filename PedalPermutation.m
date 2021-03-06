 
% Emi Fló junio 2017 @CIBPsi
%% -- Script para Permutation Cluster Analisis siguiendo el analisis realizado por (Bemis & Pylkaneen, 2011). 

% Bemis & Pylkaneen, 2011 :
%       -en ventana de tiempo 0-500 ms
%       -10.000 permutaciones
%       -Identifican cluster de tiempo adyacentes en los que haya una
%       interacción entre tarea y cantidad de palabras mediante un anova
%       2x2. Por lo menos 10 puntos de tiempo contiguos con  p<0.30.
%       - Para cada punto de tiempo de cada cluster ttest dentro de cada
%       tarea (factores una palabra vs dos palabras)
%       -Estadístico = SumaTvaluesComposicion -SumaTvaluesLista (Aumento
%       en la actividad durante la condicion dos palabras en comparación
%       con una palabra en la tarea de composición. Ninguna diff de
%       actividad entre las condiciones en la tarea lista.


%%% Este script usa las funciones:

%       -permutar_eeg
%       -clusterTiempo
%       -indicesClusters
%       -checkCluster2
%       -clusterTiempo_Espacio
%       -ttestClusterTiempoEspacio.

%%% Estas funciones se pueden usar para otros estadísticos
%%% siempre y cuando se parta de una matriz de pvalues (o lo que fuera) de ce x sf 
%%% (cantidad de electrodos x frecuencia de muestreo)

%%% Necesita array de vecinos. 

%%% Necesita data de sujeto como average, en fieldtrip salida allsubj_trig%


% Diciembre 2017 - Encuentro error en 6) y 7)- La salida de checkCluster2
% que interesa es elecSoFar. Array con matriz de cada cluster. elecSoFar
% debe ser el input de la función ttestClusterTiempoEspacio. 
% El paso 7) queda obsoleto.
%% 0) Load and prepare data

% cargo data con estructura allsubj_trig de fieldtrip. Las funciones están
% diseñadas para correr sobre la estructura de las permutaciones por lo que
% la data tendrá la misma estructura. Array de condiciones x cantidad de sujetos.

clc
clear

load ERP_25_subjects_Composition_Bl_primer_palabra_zscore101217_1410.mat
load ERP_25_subjects_ListaAdjetivo_Bl_primer_palabra_zscore101217_1341.mat
load ERP_25_subjects_ListaSustantivo_Bl_primera_palabra_zscore_101217_1407.mat

Comp2W = {};
Comp1W = {};

ListA2W = {};
ListaA1W = {};

ListN2W = {};
ListN1W = {};

ce = 64; % cantidad de electrodos
tp1 = 1; % initial time point 
tp2 = 718; % final time point
c=25; % number of subjects

for i=1:length(allsubj_trig_1WC)
    Comp2W{1,i} = allsubj_trig_2WC{1,i}.avg([1:ce],[tp1:tp2]);
    Comp1W{1,i} = allsubj_trig_1WC{1,i}.avg([1:ce],[tp1:tp2]);
    ListA2W{1,i} = allsubj_trig_2WLA{1,i}.avg([1:ce],[tp1:tp2]);
    ListA1W{1,i} = allsubj_trig_1WLA{1,i}.avg([1:ce],[tp1:tp2]);
    ListN2W{1,i} = allsubj_trig_2WLS{1,i}.avg([1:ce],[tp1:tp2]);
    ListN1W{1,i} = allsubj_trig_1WLS{1,i}.avg([1:ce],[tp1:tp2]);
end

%Real data

RealDataComp = {};
RealDataComp(1,1:c)= Comp2W(1:c);
RealDataComp(2,1:c) = Comp1W(1:c);

RealDataLA = {};
RealDataLA(1,1:c) = ListaA2W(1:c);
RealDataLA(2,1:c) = ListaA1W(1:c);

RealDataLN = {};
RealDataLN(1,1:c) = ListaS2W(1:c);
RealDataLN(2,1:c) = ListaS1W(1:c);

%% 1) Data prep

%%% Randomly permutates the condition average for each subject 
% -------------------------------------------------------------------
% Needs allsubject fieldtrip structure
% Returns a 1 x number of permutation array.

c=25; % Number of subjects
g=2000; % Number of permutations
cond = 2; % Number of conditions

perm_C = permutar_eegMatlab2017(Comp2W,Comp1W,c,g,cond);
perm_LA = permutar_eegMatlab2017(ListA2W,ListA1W,c,g,cond);
perm_LN = permutar_eegMatlab2017(ListS2W,ListS1W,c,g,cond);


%% Generate statistic distribution

save_dir = '/home/eflo/Escritorio/Analisis_Permutaciones_Composicion2.0/SalidaPermutaciones/ISIOK/';

pt = tp2-tp1+1; % Time points to evaluate
ce = 64; % Number of electrodes
c = 25;% Number of subjects

t0=tic;

load vecinos_biosemi64

for numberPerm=1:length(perm_C)
    
    tic
    
    display(numberPerm)
   
% 2)  ANOVA 2x2 for each data point

    permutation_C = perm_C{numberPerm}; 
    permutation_LN = perm_LN{numberPerm};
%     permutation_LA = perm_LA{numberPerm};
    
    FACTNAMES = {'Tarea', 'CantPalabras'};
    
    % pvalues = ce x sf 
    [ pvalues_ ] = anova2x2Permutacion_beta(c,pt,ce,FACTNAMES,...
        permutation_C,permutation_LN);

% 3) Find clusters in time
  
    threshold = 0.20; % pvalue for threshold 
    t = 5; % Number of minimum time points that have to take in part in a cluster
    umbralElectrodos = 3; %mínimo de electrodos que deben formar parte del cluster para c/punto de tiempo

    [indices_tiempo cluster_tiempo] = clusterTiempo(ce,pt,threshold,t,pvalues_);    
    [ listaElectroClusters listaClusters ] = indicesClusters(indices_tiempo, ce);
 
% 4) Encuentro clusters in space 
        
    ptable = cluster_tiempo; % 
    
    [clustSoFar,clusterMatrix, elecSoFar] = checkCluster2(ptable, vecinos,{},0,{},listaElectroClusters, listaClusters,[],0,umbralElectrodos);

% 5) Ttest for each task, for each data point and electrode taking part in a cluster
    
    alpha = 0.05;
    
    [ clusterTiempoEspacioTT_comp ] = ttestClusterTiempoEspacio(permutation_C, elecSoFar, c, ce, pt, alpha);
    
    [ clusterTiempoEspacioTT_list ] = ttestClusterTiempoEspacio(permutation_LA, elecSoFar, c, ce, pt, alpha);
    
% 6) Build statistic distribution.

    diff_C_L=zeros(1,length(clusterTiempoEspacioTT_comp));
    diff_L_C=zeros(1,length(clusterTiempoEspacioTT_comp));
    
    
    for i=1:length(clusterTiempoEspacioTT_comp)
        diff_C_L(1,i)=sum(abs(clusterTiempoEspacioTT_comp{1,i}(:)))-sum(abs(clusterTiempoEspacioTT_lista{1,i}(:)));
        diff_L_C(1,i)=sum(abs(clusterTiempoEspacioTT_lista{1,i}(:)))-sum(abs(clusterTiempoEspacioTT_comp{1,i}(:)));
    end
    
    %%% clusters : composition - list

    tclusters_all_C_L{numberPerm} = diff_C_L;
    A=find(diff_C_L>0);% only positive clusters (composition difference bigger than list difference)
    tclusters_C_L{numberPerm} = diff_C_L(A);
    
    % clusters : list - composition
       
    tclusters_all_L_C{numberPerm} = diff_L_C;
    
    B=find(diff_L_C>0);% only positive clusters (list difference bigger than composition difference)

    tclusters_L_C{numberPerm} = diff_L_C(B);
    
  toc 
  
  
end

t=toc(t0)

% load handel; sound(y,Fs) 

% save([save_dir 'filename'],'tclusters_C_L','tclusters_all_C_L','tclusters_L_C','tclusters_all_L_C');
%% Test actual data

save_dir = '/home/eflo/Escritorio/Analisis_Permutaciones_Composicion2.0/SalidaPermutaciones/ISIOK/';

pt = tp2-tp1+1; % Time points to evaluate
ce = 64; % Number of electrodes
c = 25;% Number of subjects

t0=tic;

load vecinos_biosemi64

for numberPerm=1

   
% 2)  ANOVA 2x2 for each data point

    permutacion_comp=RealDataComp;
    permutacion_lista=RealDataLN;
%     permutacion_lista=RealDataLA;

    FACTNAMES = {'Tarea', 'CantPalabras'};
    
    % pvalues = ce x sf 
    [ pvalues_ ] = anova2x2Permutacion_beta(c,pt,ce,FACTNAMES,...
        permutation_C,permutation_LN);

% 3) Find clusters in time
    
    %%% Usa funcion clusterTiempo
  
    threshold = 0.20; % pvalue for threshold 
    t = 5; % Number of minimum time points that have to take in part in a cluster
    umbralElectrodos = 3; %mínimo de electrodos que deben formar parte del cluster para c/punto de tiempo

    [indices_tiempo cluster_tiempo] = clusterTiempo(ce,pt,threshold,t,pvalues_);
    
    [ listaElectroClusters listaClusters ] = indicesClusters(indices_tiempo, ce);
 
% 4) Encuentro clusters 
        
    ptable = cluster_tiempo; % 
    
    [clustSoFar,clusterMatrix, elecSoFar] = checkCluster2(ptable, vecinos,{},0,{},listaElectroClusters, listaClusters,[],0,umbralElectrodos);

% 5) Ttest for each task, for each data point and electrode taking part in a cluster
    
    alpha = 0.05;
    
    [ clusterTiempoEspacioTT_comp ] = ttestClusterTiempoEspacio(permutation_C, elecSoFar, c, ce, pt, alpha);
    
    [ clusterTiempoEspacioTT_list ] = ttestClusterTiempoEspacio(permutation_LA, elecSoFar, c, ce, pt, alpha);

    % To obtain matrix with cluster distribution for grpahical purpose
    matriz_clusters_comp= cell(1,length(clusterTiempoEspacioTT_comp));
    matriz_clusters_list= cell(1,length(clusterTiempoEspacioTT_list));
   
    for i=1:length(clusterTiempoEspacioTT_comp)
        
        matriz_clusters_comp{1,i}= clusterTiempoEspacioTT_comp{1,i};
        matriz_clusters_lista{1,i}= clusterTiempoEspacioTT_lista{1,i};
   
    end
    
%     save([save_dir 'tclusters_Comp_Sust_2doC_TFCE_61_3e_t5_020_005'],'matriz_clusters_comp','matriz_clusters_lista');
    
% 6) Get statistic value for real data.

    diff_C_L=zeros(1,length(clusterTiempoEspacioTT_comp));
    diff_L_C=zeros(1,length(clusterTiempoEspacioTT_comp));
    
    for i=1:length(clusterTiempoEspacioTT_comp)
        diff_C_L(1,i)=sum(abs(clusterTiempoEspacioTT_comp{1,i}(:)))-sum(abs(clusterTiempoEspacioTT_lista{1,i}(:)));
        diff_L_C(1,i)=sum(abs(clusterTiempoEspacioTT_lista{1,i}(:)))-sum(abs(clusterTiempoEspacioTT_comp{1,i}(:)));
    end
 
  % cluster : composition - list 
    
    tclusters_real_data_all_C_L = diff_C_L; 
    A=find(diff_C_L>0); 
    tclusters_real_data_C_L = diff_C_L(A);

% % cluster : list - composition
    tclusters_real_data_all_L_C = diff_L_C;
    B=find(diff_L_C>0);
    tclusters_real_data_L_C = diff_L_C(B);
end

 
% save ([save_dir 'RealData_Comp_Sust_BPP_2doC_TFCE_61_3e_t5_020_005'],'tclusters_real_data_C_L','tclusters_real_data_all_C_L','tclusters_real_data_L_C','tclusters_real_data_all_L_C')

%% Test for significant clusters

%Comp vs List

distribucion_clusters = [];

h=1;
for i=1:length(tclusters_C_L)
    for j=1:length(tclusters_C_L{i})
        distribucion_clusters(1,h)=tclusters_C_L{i}(j);
        h=h+1;
    end
end
toc

hist(distribucion_clusters)

Media = mean(distribucion_clusters);
St = std(distribucion_clusters);

cluster_significativo=[];

j=0;
for i=1:length(tclusters_real_data_C_L)
    [h,p,ci,zval]=ztest(tclusters_real_data_C_L(i),Media,St,'tail','right')
    if h==1
        cluster_significativo(j+1,1:2) = [i p];
        j=j+1;
    end
end

%% Test for significant clusters

%%List vs composition


tic
distribucion_clusters = [];

h=1;
for i=1:length(tclusters_L_C)
    for j=1:length(tclusters_L_C{i})
        distribucion_clusters(1,h)=tclusters_L_C{i}(j);
        h=h+1;
    end
end
toc

hist(distribucion_clusters)

Media = mean(distribucion_clusters);
St = std(distribucion_clusters);

cluster_significativo=[];

j=0;
for i=1:length(tclusters_real_data_L_C)
    [h,p,ci,zval]=ztest(tclusters_real_data_L_C(i),Media,St,'tail','right')
    if h==1
        cluster_significativo(j+1,1:2) = [i p];
        j=j+1;
    end
end


