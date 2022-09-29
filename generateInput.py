"""
Generate input for model.
"""

import json
import numpy as np



#user inputs:
jsonFileName = "build/inputs.json"

antennaFreq = 51.0e6 #Hz
elemCharge = 1.602176634e-19

#convert to python dict
dataDict = {
        "resolution" : 300,
        "R_west" : 2.3,
        "R_east" : 4.0,
        "R_ant" : 3.9,
        "nToroidal" : 26,
        "R0" : 3.0,
        "B0" : 3.45,
        "omega" : 2.0*np.pi*antennaFreq,
        "J0R" : 0.0,
        "J0Phi" : 0.0,
        "J0Z" : 1.0,
        "saveMatrix" : False,
        "solveMatrix" : True,
        "HarmonicMin" : -3,
        "HarmonicMax" : 3,
        "electronDensity" : 3.0e19,
        "minElecDensity" : 3.0e17,
        "hydrogenConcentration" : 0.05,
        "electronTemp" : 3.0e3*elemCharge,
        "minTemp" : 1.0e2*elemCharge,
        "angularResolution" : 101,
        "method" : 2}

#convert to JSON
jsonDict = json.dumps(dataDict, indent=2, sort_keys=True)

#write to disk
with open(jsonFileName, "w") as text_file:
    text_file.write(jsonDict)

#finished  
print("json file created")
