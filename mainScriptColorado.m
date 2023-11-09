
clear
close all

fileName = "colorado";
alphas = [ 0.1, 0.25, 0.5]; 

for alpha = alphas

    filePath = "./data/" + fileName;
    data = load( filePath );
    
    [ data, settings] = preprocessColorado( data );
    settings.concentrationParameter =  alpha * settings.numberOfExperts;

    resultObject = smc2gpmoe( data, settings);
    saveResults( resultObject, fileName);
end