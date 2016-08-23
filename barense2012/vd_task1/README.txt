27feb2016
ps

noPreTrain_new_par -> an attempt to make a much much simpler model, one without
peak and total activations, also without a complex activation function


**************************************************************

1jan2016
ps

any folder with 'tripleDims' in name refers to how the input dimensions tripled, as per
Cowell et al. 2006's actual stimuli (rather than caudal layer receiving 2 dims as input, and PRC 8,
they received 6 and 24, respectively. Though, the dims were just triplicaitons of what was 
already being input (i.e., [.35,.95] -> [.35,.35,.35,.95,.95,.95]))


**************************************************************
6jan2016

ps

NOTE: all previous earlyHalt files were from files incorrectly coded! 
The initial weights were being assigned incorrectly


21dec2015

ps

earlyHalt: prevents excessive updating of weights after intial selectivity (sevond stim) says 'novel!
     NOTE: weird thing happening in which peak activations are currently decreasing. not sure why...

nosieOnFamil: more akin to classic signal detection theory in which there is a distribution of foil and target

December 14, 2015

ps

vd_task2_updated is the version of the experiments 3/4 that were updated with the modifications
made to vd_task1 (in which famil Diff was allowed to go below 0)

October 31 2015

ps

FOLDER INFORMATION

vd_task1 -> models Barnese et. al expt 1, in which performance drops in
            the second half of trials for the lesioned network, HA only
        CURRENTLy -> parameters work, but only when 5 networks are modelled at once
			->->error bars are too small when running more rats

vd_task2 -> to model Barnese et. al expt2, in which performance is fine for the 
		first third of trials (uses pictures as well as abstract shapes, LA, 16 stims)
		drops in the second third (HA, abstract shapes only, 6 stims)
		rises again in the final third (HA again, 16 stims)

