import ROOT
import root_numpy
import pandas as pd

file    = ROOT.TFile('NuMI_Run1_Overlay_Output.root')
tree    = file.Get('UBXSec/pottree')

arr     = root_numpy.tree2array(tree, branches = ['pot'])
df      = pd.DataFrame(arr)

totpot     = 0

for i in range(df.shape[0]):

    totpot     = totpot + df['pot'][i]
    

print "The amount of total POT = %d." %( totpot )

