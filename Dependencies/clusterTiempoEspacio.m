function [ clusterTiempoEspacio ] = clusterTiempo_Espacio( clustSoFar )


%%% Funcion para generar, a partir de array clustSoFar obtenido de la
%%% función checkCluster2, una matriz cexsf para cada cluster. Los unos
%%% representan puntos de tiempo/electrodo que forman parte de ese cluster


  clusterTiempoEspacio = {}; %array que contiene una matriz ce x sf por cluster
    
    
    for i=1:length(clustSoFar)
        cluster = zeros(ce,sf);
        for j=1:length(clustSoFar{i})
            if clustSoFar{1,i}(j,1)~=0
                cluster(clustSoFar{1,i}(j,1),clustSoFar{1,i}(j,2):clustSoFar{1,i}(j,3)) = 1;
            end
        end
    clusterTiempoEspacio{1,i}= cluster;
    end
    
end

