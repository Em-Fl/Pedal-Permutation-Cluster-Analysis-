%%%  Función para recorrer los electrodos buscando clusters multielectrodo de p-valores
%%%  -------------------------------------------------------------------
%   
%  Función recursiva: algunos inputs deben ser fijados vacíos.
%
%  entradas: electrodo: electrodo a chequear (fijar vacio primera vez)
%            ptable: salida de estadístico punto a punto
%            indiceCluster1E: indices de cluster de 1 electrodo
%            chanList: lista de índices electrodos (ojo recursivo)
%            vecinos: lista de indices de electrodos vecinos
%            clustSoFar, a b: vacío la primera vez
%            esteCluster: vacío la primera vez
function checkElectrode(esteElectrodo, ptable, indiceCluster1E, chanList, vecinos, clustSoFar, esteCluster) 
    
    %inicializar variables para la primera llamada a la función
    if isempty(esteElectrodo)
        [nElec,b]=size(ptable);
        a=1;
        esteElectrodo=1;
        electrodos=1:nElec;
        clustSoFar=zeros(nElec,b);
    end
    
    % para cada cluster
    for i=1:length(indiceCluster1E(esteElectrodo))
        [a,b]=indiceCluster1E(esteElectrodo,i); % ver implementación de esto
        %chequear si este cluster ya forma parte de un cluster
        if clustSoFar[esteElectrodo,a:b]~=0
            
        end
        %multielectrodo
        
    end
    
    
    return(electrodos, a,b)
end