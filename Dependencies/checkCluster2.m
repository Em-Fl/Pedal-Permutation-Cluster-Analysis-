%%%  Función para recorrer los clusters buscando clusters multielectrodo de p-valores
%%%  Álvaro Cabana Junio 2017 - CIBPsi
%%%  -------------------------------------------------------------------
%   
%  Función recursiva: algunos inputs deben ser fijados vacíos.
%  
%
%  entradas: ptable: salida de estadístico punto a punto [1s y 0s]
%            vecinos: cell array, para cada electrodo, vector con indices de electrodos vecinos
%            clustSoFar: vacio la primera vez - 
%            elecSoFar: vacio la primera vez - electrodos en este cluster, hasta ahora
%            esteCluster: vacío la primera vez - 
%            listaElectroClusters: matriz nClusters x 3 (E, a ,b)
%            listaClusters: cell array de electrodos, cada cell lista de
%               indices de principio y fin de cluster (a,b)
%            numCluster: uso interno, poner 0!!
%            umbralElectrodos - número mínimo de electrodos que tiene que
%            tener el cluster TODO el TIEMPO (hard threshold)
%      esteCluster tiene formato "electrodo . a,b"
%
%
%   salida final:  clusSoFar: dummy output
%                  clusterMatrix: dummy output
%                   elecSoFar: [1 nClusters] cell array con clusterMatrix para cada
%                       cluster (nElec, t)

function [clustSoFar,clusterMatrix, elecSoFar]=checkCluster2(ptable, vecinos, clustSoFar, elecSoFar, esteCluster, listaElectroClusters, listaClusters, clusterMatrix,numCluster,umbralElectrodos) 
    [nElec,t]=size(ptable);
    % Primera corrida, esteCluster vacio
    if isempty(esteCluster)
        numCluster=0;
        %inicializo matriz de clusters
        clusterMatrix=zeros(nElec,t);

    
        % aquí voy a juntar TODOS los macroClusters encontrados:
        macroCluster={};
        
        % aquí voy a juntar todos los clusters del PRESENTE macroCluster
        %------------

        %chequeo que haya algún
        
        % De la lista de clusters, agarro el primer cluster 
        
        %chequeo si algún electrodo tiene cluster por las dudas
        for i=1:length(listaElectroClusters(:,1))
            esteCluster=listaElectroClusters(i,:);
            electrodo=esteCluster(1); % agarro el primer cluster del primer electrodo que tenga cluster
            a=esteCluster(2);
            b=esteCluster(3);
            
            %chequeo si está en algún cluster de algún lado so far 
            if clusterMatrix(electrodo,a)==0
                numCluster=numCluster+1;
                clustSoFar=zeros(10,3);
                elecSoFar=0;
                [clustSoFar,clusterMatrix, elecSoFar]=checkCluster2(ptable, vecinos, clustSoFar, elecSoFar, esteCluster, listaElectroClusters, listaClusters, clusterMatrix,numCluster,umbralElectrodos);
                macroCluster = [macroCluster  {clustSoFar}];
                
            end
        end
        
        % para cada cluster, en clustSoFar, obtener t minde columna 2 y t
        % max de columna 3
        % ir a MAtrixcluster, meter suma de columnas entre esos t, y
        % chequear que la suma sea igual o mayor a el parámetro
        % numeroVecinos
        % filtrar si no
        
        clustSoFar=macroCluster;
        % vamos a quedarme con los clusters que tienen menos de
        % umbralElectrodos
        quedanClusters=true(1,numCluster);
        arrayMatrices={};
        for i=1:numCluster
            clusMxi=clusterMatrix==i;
            sumaMx=sum(clusMxi);
            cluster=sumaMx>0;
            %suma de Puntos por columna (electrodos por tiempo)
            llegoAUmbral=sumaMx(cluster)>=umbralElectrodos;
            noGo=false;
            %inicio de Cluster
            indicesCluster=find(cluster);
            inicio=indicesCluster(1);
            % si no hay punto de tiempo alguno con igual o mayor numero de
            % electrodos a umbral, chim pum fuera
            if sum(llegoAUmbral)==0
                noGo=true;
            else
                % recortar las puntas (onsets y offsets de cada cluster)
                posiciones=find(llegoAUmbral);
                onset=posiciones(1);
                offset=posiciones(end);
                % no puede haber ceros en llegoAUmbral entre el primer y el
                % ultimo 1
                if sum(~llegoAUmbral(onset:offset))
                    noGo=true;
                else
                    clusMxi(:,1:(inicio+onset-1))=0;
                    clusMxi(:,(inicio+offset):end)=0;
                end
                % onset    
            end
            
            if noGo
                quedanClusters(i)=0;
                clusterMatrix(clusMxi)=0;
            else
                thisClusterMatrix=zeros(nElec,t);
                thisClusterMatrix(clusMxi)=1;
                arrayMatrices=[arrayMatrices thisClusterMatrix];
            end
        end
        clustSoFar=clustSoFar(quedanClusters);
        elecSoFar=arrayMatrices;
 
    else
        a=esteCluster(2);
        b=esteCluster(3);
        esteElectrodo=esteCluster(1);

        %agrego este cluster a una lista de clusters multiElectrodo....
        %tic
        %clustSoFar=[clustSoFar esteCluster]; %ineficiente pero robusto (?!
        
        % si hay 10, 20, 30 electrodos... etc. en el cluster, agrego 10 màs
        % a la matriz de clustSoFar
        if elecSoFar>0 && (mod(elecSoFar,10)==0)
            clustSoFar = [clustSoFar ; zeros(10,3)];
        end
        % un electrodo mas al cluster
        elecSoFar=elecSoFar+1;
        clustSoFar(elecSoFar,:)=esteCluster;
        %toc
        % agregamos el cluster a la matriz de clusterMatrix (para chequear
        % rapido)
        clusterMatrix(esteElectrodo,a:b)=numCluster;  %quick and dirty: ver chequeo cluster


        % ir a todos los electrodos vecinos de este electrodo, y mirar ahi si
        % solapan sus clusters
        estosVecinos=vecinos{esteElectrodo};
        % para todos los vecinos
        for i=1:length(estosVecinos)

            otroElectrodo=estosVecinos(i);
            losClusters=listaClusters{otroElectrodo};
            if ~isempty(losClusters)
                for j=1:length(losClusters)
                    otroCluster=losClusters{j};

                    a2=otroCluster(1);
                    b2=otroCluster(2);

                     %chequear que no sea un cluster previamente visitado
                     % si no lo visité, entonces sigo adelante
                    if clusterMatrix(otroElectrodo,a2)==0               
                        % CHEQUEOOOOO   QUE  se solapen!
                        if (( a2<b && a<b2 ) || (a<b2 && a2<b ) )
                            [clustSoFar,clusterMatrix, elecSoFar]=checkCluster2(ptable, vecinos, clustSoFar, elecSoFar,[otroElectrodo otroCluster], listaElectroClusters, listaClusters, clusterMatrix,numCluster,umbralElectrodos);
                        end
                    end
                end
            end
        end       
    end

end

