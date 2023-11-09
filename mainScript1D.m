
clear
close all

fileNames = ["stationary", "nonstationary", "motorcycle"];
alphas = [ 0.1, 1, 3.5];

nFiles = length( fileNames );

for fileName = fileNames
    for alpha = alphas

        filePath = "./data/" + fileName;
        data = load( filePath );
        
        [ data, settings] = preprocess1D( data );

        settings.showDebugFigure = true;
        settings.concentrationParameter =  alpha;

        resultObject = smc2gpmoe( data, settings);
        saveResults( resultObject, fileName);
    end
end

