
function [indices_tiempo cluster_tiempo] = clusterTiempo(ce,sf,threshold,t,pvalues_)

%%% Funcion que encuentra los puntos de tiempo consecutivos que conforman
%%% un cluster en cada electrodo, de acuerdo a parámetros establecidos.
%%% Hay que definir valor de p umbral y cantidad mínima de puntos
%%% consecutivos con ese valor de p para conformar un cluster temporal
%%%
%%% La salida es 
%%% - indices_tiempo: array de 1 x 64 que contiene para cada 
%%% electrodo un cell array con vectores que corresponde con los indices 
%%% de tiempo de los clusters
%%% - cluster_tiempo: matriz de ce x sf con unos marcando los puntos de 
%%% tiempo/espacio donde el pvalue es menor al umbral establecido.

% *uso funciones de Image Processing toolbox

% Argumentos que precisa la funcion

% sf frecuencia de muestreo 
% ce cantidad de electrodos
% threshold valor p umbral de significancia
% t cantidad de puntos de tiempo mínimos que deben conformar un cluster
% pvalues matriz de ce x sf con los pvalues obtenidos del estadístico

largoSpan = zeros(ce, sf);
spanLocs = zeros(ce, sf);
spanLength = {};
area = {};
goodSpans = {};
allInSpans = {};
indices_tiempo = {}; %array de 1 x 64 que contiene para cada electrodo un cell array con vectores que corresponde con 
%los indices de tiempo de los clusters
cluster_tiempo = zeros(ce,sf);


for i=1:ce
    belowThreshold = (pvalues_(i,:) < threshold);
    spanLocs(i,:) = bwlabel(belowThreshold);   %identify contiguous ones
    spanLength{i} = regionprops(spanLocs(i,:), 'area');
    area{i} = [spanLength{i}.Area]; %length of each span
    goodSpans{i} = find(area{i}>=t); %get only spans of t+ time points
    allInSpans{i} = find(ismember(spanLocs(i,:), goodSpans{i}));  %indices of these spans uno a continuaci�n del otro
    for v=1:length(goodSpans{1,i})
        indices_tiempo{i}{1,v} = find(spanLocs(i,:) == goodSpans{i}(v));
    end
    if i==64 && isempty(goodSpans{1,i})
        indices_tiempo{64}=[];
    end
end


for l=1:ce
    for m=1:length(allInSpans{1,l})
        cluster_tiempo(l,(allInSpans{1,l}(m))) = 1;
    end
end


