function [ pvalues_ ] = anova2x2Permutacion(c,sf,ce,FACTNAMES,permutacion_comp,permutacion_adj )

%%% Utiliza la funcion rm_anova2 descargada de mathworks (https://www.mathworks.com/matlabcentral/fileexchange/6874-two-way-repeated-measures-anova?s_tid=srchtitle) 
%%% La salida es una matriz de dimensiones cantidad de electrodos x
%%% frecuencia de muestreo

%%% Los argumentos de la función son:
    % - c : cantidad de sujetos
    % - sf: frecuencia muestreo
    % - ce: cantidad de electrodos
    % - FACTNAMES : Nombre de los factores
    % - permutaciones de cada tarea con la estructura allsubj salida de
    % fieldtrip

pvalues_ = zeros(ce, sf);

for t=1:sf
    
    for e=1:ce
        
        et = []; 
        F1 = []; %Factor1 = TAREA (1.Composicion, 2.Lista)
        F2 = []; %Factor2 = cantidad de palabras (1.OW 2.TW)
        for i=1:c
            et((length(et)+1),1)= permutacion_comp{1,i}(e,t);
            F1((length(F1)+1),1) = 1;
            F2((length(F2)+1),1) = 2;
        end
        for i=1:c
            et((length(et)+1),1)= permutacion_comp{2,i}(e,t);
            F1((length(F1)+1),1) = 1;
            F2((length(F2)+1),1) = 1;
        end
        for i=1:c
            et((length(et)+1),1)= permutacion_adj{1,i}(e,t);
            F1((length(F1)+1),1) = 2;
            F2((length(F2)+1),1) = 2;
        end
        for i=1:c
            et((length(et)+1),1)= permutacion_adj{2,i}(e,t);
            F1((length(F1)+1),1) = 2;
            F2((length(F2)+1),1) = 1;
        end
        
        %hago el vector columna que identifica cada sujeto
        s = [];
        B = 1:c;
        s =[B B B B];
        S=s';

        Y = et;

        %anova de medidas repetidas
        salida = rm_anova2(Y,S,F1,F2,FACTNAMES);
        p = salida{4,6}; %pvalue de la interaccion tarea x cant_clust_tiempo de palabras

        %guardo los valores de p para cada punto de tiempo para cada
        %electrodo
        pvalues_(e,t) = p;
    end
end



