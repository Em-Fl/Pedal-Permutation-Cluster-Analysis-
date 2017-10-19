function [ pvalues_ ] = anova2x2Permutacion_matlab2017(c,sf,ce,permutacion_comp,permutacion_adj )

%%% Este comentario en stackoverflow me salvó la vida para lograr la salida
%%% deseada con la función fitrm y ranova https://www.mathworks.com/matlabcentral/answers/140799-3-way-repeated-measures-anova-pairwise-comparisons-using-multcompare

%%% Los argumentos de la funciÃ³n son:
    % - c : cantidad de sujetos
    % - sf: frecuencia muestreo
    % - ce: cantidad de electrodos
    % - FACTNAMES : Nombre de los factores
    % - permutaciones de cada tarea con la estructura allsubj salida de
    % fieldtrip. Lo armé para que que en la fila 1 de los arrays esté la
    % condicion dos palabras, y en la fila 2 la condicion 1 palabra. Hay un
    % array distinto para cada condicion de tarea.

pvalues_ = zeros(ce, sf);

for t=1:sf
    
    for e=1:ce
        
        comp2w=zeros(23,1);
        comp1w=zeros(23,1);
        adj2w=zeros(23,1);
        adj1w=zeros(23,1);
        
        for i=1:c
            
            comp2w(i,1) = permutacion_comp{1,i}(e,t);
            comp1w(i,1) = permutacion_comp{2,i}(e,t);
      
            adj2w(i,1) = permutacion_adj{1,i}(e,t);
            adj1w(i,1) = permutacion_adj{2,i}(e,t);
            
        end
        
        %hago el vector columna que identifica cada sujeto
        s=[1:c];
        S=s';

        %Genero tabla con columna que identifica los sujetos y los valores
        %de voltaje para cada condicion
        datos = table(S,comp2w,comp1w,adj2w,adj1w,'VariableNames',{'sujeto','comp2w','comp1w','adj2w','adj1w'});
       
        factores = [1 2;1 1;2 2;2 1];
        factoress=table(factores(:,1),factores(:,2),'VariableNames',{'Tarea','cantidadPalabras'});
        
        Factores =factoress;
        %convierto a categoricos
        Factores.Tarea=categorical(Factores.Tarea);
        Factores.cantidadPalabras=categorical(Factores.cantidadPalabras);
        Factores.interaccion=Factores.Tarea .* Factores.cantidadPalabras; %factor de la interacción
        

        rm = fitrm(datos,'comp2w-adj1w~1','WithinDesign',Factores);
        ranovatbl = ranova(rm);
        
        [ranovatbl] = ranova(rm,'WithinModel','Tarea*cantidadPalabras');
        
        pvalues_(e,t)= table2array(ranovatbl(7,5));

    end
end
