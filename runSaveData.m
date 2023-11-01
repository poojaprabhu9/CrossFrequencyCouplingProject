% runSaveData

folderSourceString = 'N:\Projects\Pooja_CrossFrequencyCouplingProject\';
monkeyName = 'alpaH'; expDate = '210817'; protocolName = 'GRF_002';
removeEvokedResponse=1; tapers = [1 1];

saveData(monkeyName,expDate,protocolName,folderSourceString,removeEvokedResponse,tapers);
