function [ pvalues_] = ttestPareadoParaValidarConFieldtrip(permutacion, c, ce, sf, alpha)
%Función que realiza un ttest pareado sobre la estructura permutado.
%Permutado es un cellarry de 1xcantidad de permutaciones. Dentro de cada cell
%array se encuentran un cellarray de cantidad de condiciones x cantidad de
%sujetos.

% c cant de sujetos
% ce cant de electrodos
% sf frecuencia de muestreo

pvalues_ = zeros(ce, sf);
% tvalues = zeros(ce,sf);

%En cada vuelta et es un vector de voltaje para un punto de tiempo para un
%electrodo de largo cant de sujetos, cada valor es el voltaje para ese 
%punto de tiempo para cada sujeto

for e=1:ce
    for t=1:sf
        et_2WC = zeros(c,1);
        et_1WC = zeros(c,1); 
        
        for i=1:c
            
            et_2WC(i,1)= permutacion{1,i}(e,t);
            et_1WC(i,1)= permutacion{2,i}(e,t);
            
        end
            
            %ttest para cada punto de tiempo entre condicion una palabra vs
            %dos palabras
            
            [~,p,~,stats] = ttest(et_2WC,et_1WC,alpha,'both');
        
            pvalues_(e,t) = p;
            %tvalues(e,t) = stats.tstat;
            
    end
    
end


