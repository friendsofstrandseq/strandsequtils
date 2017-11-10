In comparison to git version:
added readbamsAndNBfit.R
added SVcalling_cmd.R
added runScript.R
replaced temp with outputDir in all writing outout lines (in readbamsAndNBfit.R in comparison to call.R)-- should also be fixed in call.R
added newSVcalling.R, newSVcalling_cmd.R, newhapStatusCellProbTable.R, getSisterHaplotype.R: IN these three files I fixed the bug in computing aggregate probabilities, so they should be replaced with the old ones
added runClassificationOnSimulatedData.R --> run the whole thing on Sascha's simulated data

