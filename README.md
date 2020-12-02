# uBooNE-NuMI-NuMuCCInclusive-From-Target-Contained-Analysis
This is the analysis code used for my Ph.D. thesis measurement, the differential cross section of contained muon neutrino events from the NuMI beam target with MicroBooNE.

This analysis involves a set of box cuts designed to separate out background from signal.  No machine learning framework is used in this analysis's selection.

This repository consists of three files:

1. UBXSec_module.cc - This file performs part of the selection and fills the output "passing_events_tree" with the information necessary to complete the selection.  The 'Truth_Tree' contains truth information of the signal events input into the selection which used to calculate efficiency as a function of muon truth kinetic energy and angle.

2. Making_Cut_By_Cut_Plots.C - This file performs the second part of the selection and plots a number of the output quantities of events passing the full selection.  It is written to easily turn parts of the selection off to see what effect(s) that has on agreement between data and prediction.

3. counting_mc_pot.py - This file counts the amount of Protons On Target (POT), the metric used to normalized the samples, in the simulation.
