function [ listaElectroClusters listaClusters ] = indicesClusters(indices_tiempo, ce)

 %%% Funcion para generar dos cellarrays necesarios para la funcion checkCluster2
 
 %%% 1) listaElectroClusters matriz nClusters x 3 (E, a,b)
    
    % E=Electrodo donde está el cluster
    % a=indice de tiempo donde comienza cluster
    % b=indice de tiempo donde termina cluster
    
 %%% 2)  listaClusters 
    %%% cell array de electrodos, cada cell lista de indices de principio y fin de cluster (a,b).

 %%% Argumentos que requiere la función
      % -  indices_tiempo - array 1x 64 que contiene para cada electrodo un cell array con vectores con indices de tiempo de los clusters
      % - ce = cantidad de electrodos
      
    listaElectroClusters = [];
    
    l=0;
    
    for i=1:length(indices_tiempo)
        if ~isempty(indices_tiempo{1,i}) % solo si hay cluster en el electrodo
            for j=1:length(indices_tiempo{1,i})
                l=l+1;
                %indices_tiempo{1,i}{1,j}
                listaElectroClusters (l,:) = [i indices_tiempo{1,i}{1,j}(1) indices_tiempo{1,i}{1,j}(end)];
            end
        end
    end
    
    listaClusters = {};
    
    %Requiere indices_tiempo - array 1x 64 que contiene para cada electrodo un cell array con vectores con
    %los indices de tiempo de los clusters
    
    for i=1:ce
        if  ~isempty(indices_tiempo{1,i})
            
            for j=1:length(indices_tiempo{1,i})
                
                listaClusters{1,i}{j}=[indices_tiempo{1,i}{1,j}(1),indices_tiempo{1,i}{1,j}(end)];
                
            end
        else
            listaClusters{1,i} = {};
        end
    end

end

