function [ pvalues_ ] = anova2x2Permutacion_beta(c,sf,ce,FACTNAMES,permutacion_comp,permutacion_adj )

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

for e=1:ce  %electrodo 1
    
    et = []; 
    F1 = []; %Factor1 = TAREA (1.Composicion, 2.Lista)
    F2 = [];
    S=[];
    nCond=4;

    for s=1:c
        i=nCond*(s-1);
        S(i+1)=s;
        et(i+1,:)= permutacion_comp{1,s}(e,:);
        F1(i+1) = 1;
        F2(i+1) = 2;
        et(i+2,:)= permutacion_comp{2,s}(e,:);
        F1(i+2) = 1;
        F2(i+2) = 1;
        S(i+2)=s;
        S(i+3)=s;
        et(i+3,:)= permutacion_adj{1,s}(e,:);
        F1(i+3) = 2;
        F2(i+3) = 2;
        et(i+4,:)= permutacion_adj{2,s}(e,:);
        F1(i+4) = 2;
        F2(i+4) = 1;
        S(i+4)=s;
    end

        Y = et;

        %anova de medidas repetidas
        salida = rm_anova2_beta(Y,S,F1,F2,FACTNAMES);
        p = salida{4,6}; %pvalue de la interaccion tarea x cant_clust_tiempo de palabras

        %guardo los valores de p para cada punto de tiempo para cada
        %electrodo
        pvalues_(e,:) = p;
end




