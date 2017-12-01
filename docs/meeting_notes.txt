********************************************
20171114 Meeting notes (AJB, HH, KC, HB, LM)
********************************************
PURPOSE: determine the error probability of filling kspace with an anomalous beat during ECG-gated 2D PCMRI

NOTES/ACTION ITEMS:

- peform lit review (20171201 update: ongoing)
	o HH: center kspace issues (start with Markl 2001 and do a forward lit search)
	o AJB: variability in flow in vivo (pulmonary flow most likely)
	
- simulate time-resolved data and kspace filling probabilities
	x AJB: to send out initial 'feasibility' code (put on SVN?  20171201 update- now on github)
		x done, see repo here https://github.com/alexbarks/beat_variability
	o HB/LM/HH: time domain, vary kspace filling - probability of significant error? (>10%)
		* compute error bounds (what is 95% confidence interval?)
		* 20171201 update: ongoing, LM to share progress on server here:
		  \\10.254.136.37\data_imaging\cv_mri\beat_variability\_progress_20171201\
	
- in vivo realtime measurements (to get a feel for variability)
	o AJB/HH: pulmonary/aortic real time measurements in volunteers?
	o KC: aquire some real time data on voluteers (time permitting)

********************************************
20171201 Meeting notes (AJB, HB, HH, LM)
********************************************

UPDATES/NOTES/ACTION ITEMS

- Continue lit review (all)... see 2017114 for additional notes
	o HH (and others) to put papers on server here:
	  \\10.254.136.37\data_imaging\cv_mri\beat_variability\_literature\
	  
- continue efforts from 20171114 to simulate probability of significant error
	o HB/HH/LM: cycle through different starting points and linearly filling kspace
	
- New questions
	o AJB: is 2D PCMRI typically breathhold or respiratory gated (what is acquisition time)?
		o AJB typical matrix of the acquisition?
	
FUTURE PROJECT IDEAS
20171114
- 4D flow, multiple central kspace lines (beat beat) (KC)
- collect inverse RECAR? at the same time as +RECAR  (KC)
- detect abnormal flow during aquisition (somehow automatically track the center of kspace and flow position)?
