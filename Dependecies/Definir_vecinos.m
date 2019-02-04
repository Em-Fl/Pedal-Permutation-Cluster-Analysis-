%% Script para definir vecinos para hacer Cluster Permutation Analysis con
%% script pedal_permutation de Emi Fló y Álvaro Cabana junio 2017 @CIBPsi

% Necesito el template de vecinos de biosemi y algún dato del que obtener
% la información asociada a el nombre y posición del electrodo en la
% matriz

load biosemi64_neighb

sf=512; %frecuencia de muestreo
ce=64; %cant_clust_tiempo de electrodos

vecinos = {}; %Array que contiene en orden de acuerdo al layout de biosemi 
%los electrodos. El electrodo 1 (posicion 1 en array, corresponde con Fpz) del array tiene como
%vecinos al 2,3,33,37

electrodos = permutado{1,1}{1,1}.elec.label;
for e=1:ce
    for h=1:length(neighbours(1,e).neighblabel)
        vecinos{e,1}(h) = strmatch(neighbours(1,e).neighblabel(h), electrodos, 'exact');
    end
    %vecinos{e,1} = [e, vecinos{e,1}] %por si quisiera tener el electrodo y
    %sus vecinos 
end

save('vecinos_biosemi64','vecinos')
