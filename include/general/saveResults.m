function [] = saveResults( smc2Object, fileName)

    settings = smc2Object.settings;
    folderName = fileName + "_" + num2str( settings.concentrationParameter );

    timeStamp = datetime('now', 'Format', 'yyyy-MM-dd-HH-mm-SS' );
    timeStamp = string( timeStamp );

    saveFolder = "./results/" + folderName + "_" + num2str(timeStamp);
    saveFolder = strrep( saveFolder, '.mat', '');
    mkdir( saveFolder );

    fullSaveFileName = saveFolder + "/resultObject.mat";
    save( fullSaveFileName, '-struct', 'smc2Object', '-v7.3');
end

