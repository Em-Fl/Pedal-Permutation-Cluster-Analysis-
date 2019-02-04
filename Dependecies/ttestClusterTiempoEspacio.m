function [ clusterTiempoEspacioTT ] = ttestClusterTiempoEspacio( permutacion, clusterTiempoEspacio, c, ce, sf, alpha)

%%% Funcion que hace un ttest para los cluster tempoespaciales encontrados.
%%% Para cada punto de tiempo/electrodo que forma parte de un cluster
%%% obtengo un valor de t. La salida es un cell array con 1 x cantidad de 
%%% clusters, cada cell contiene una  matriz de ce x sf con los valores de 
%%% t que forman parte de cada cluster.

%%% Argumentos que requiere la funcion:
    % - clusterTiempoEspacio: cell array obtenido de funcion c

clusterTiempoEspacioTT = {};
    
    for i=1:length(clusterTiempoEspacio)
        
        [elec tiempo] = find(clusterTiempoEspacio{1,i}); %encuentro para cada cluster los puntos de tiempo distintos de cero.
        %la salida son dos
        %vectores (electrodo y
        %tiempo) con los índices
        %de los puntos que forman
        %parte del cluster
      
        tvalues_ = zeros(ce,sf);
        
        for j=1:length(elec)
            
            et_2WC = zeros(c,1);
            et_1WC = zeros(c,1);

            
            for h=1:c %loopea entre sujetos para levantar los valores de voltaje para ese punto de tiempo
                %el ttest lo hago a partir de un vector de voltajes promedio de
                %cada sujeto para cada una de las condiciones (et_2wc y
                %et_1WC)
                
                et_2WC(h,1)= permutacion{1,h}(elec(j),tiempo(j));
                et_1WC(h,1)= permutacion{2,h}(elec(j),tiempo(j));
                
            end
            
            %ttest dentro de cada tarea, condicion una o dos palabras
            
            [~,~,~,statsC] = ttest(et_2WC,et_1WC,alpha,'both'); %demora por stat!
            tvalues_(elec(j),tiempo(j)) = statsC.tstat;
            
            clusterTiempoEspacioTT{1,i} = tvalues_;
            
        end
    end
    
end

