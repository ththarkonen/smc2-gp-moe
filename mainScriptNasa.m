
clear
close all

fileName = "nasa";
alphas = [ 0.1, 0.25, 0.5]; 

for alpha = alphas

    filePath = "./data/" + fileName;
    data = load( filePath );
    
    [ data, settings] = preprocessNasa( data );
    settings.concentrationParameter =  alpha * settings.numberOfExperts;

    resultObject = smc2gpmoe( data, settings);
    saveResults( resultObject, fileName);
end