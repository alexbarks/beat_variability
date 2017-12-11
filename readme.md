
## 20171114 Meeting notes (AJB, HH, KC, HB, LM)

**PURPOSE**: determine the error or probability of filling central kspace with an anomalous beat during ECG-gated 2D PCMRI

**NOTES/ACTION ITEMS:**

- [ ] peform lit review (20171201 update: ongoing)
	- [ ] @curatv: center kspace issues (start with Markl 2001 and do a forward lit search)
	- [ ] @alexbarks: variability in flow in vivo (pulmonary flow most likely)
	
- [ ] simulate time-resolved data and kspace filling probabilities
	- [x] ~@alexbarks: to send out initial 'feasibility' code (put on SVN?  20171201 update-now on github here)~
	- [x] ~HB/LM/HH: time domain, vary kspace filling - probability of significant error? (>10%)~
	- [ ] compute error bounds (what are confidence intervals? Monte Carlo simulation?)
	
- [ ] in vivo realtime measurements (to get a feel for variability)
	- [ ] AJB/HH: pulmonary/aortic real time measurements in volunteers?
	- [ ] KC: aquire some real time data on voluteers (time permitting)
	- [ ] AJB: does Lon or Ning have data already?

## 20171201 Meeting notes (AJB, HB, HH, LM)

**UPDATES/NOTES/ACTION ITEMS**

- [ ] All: sign up for github account and request to be a collaborator (or have alex invite you)
- [ ] All: upload your latest code and let's figure out how to delegate tasks
- [ ] All: Continue lit review... see 2017114 for additional notes
	- [ ] HH (and others) to put papers on server here:
	  \\10.254.136.37\data_imaging\cv_mri\beat_variability\_literature\

- [x] AJB: is 2D PCMRI typically breathhold or respiratory gated (also what is acquisition time/matrix)?
	* if respiratory-gated, then usually ~4 lines per segment
	  
- [ ] ongoing efforts from 20171114 action item to simulate probability of significant error
	- [ ] HB/HH/LM: cycle through different starting points and linearly fill kspace?
	- [ ] LM/HH/HB: share progress on server here:	
	  \\10.254.136.37\data_imaging\cv_mri\beat_variability\_progress_20171201\
	
## FUTURE PROJECT IDEAS
- 20171114, 4D flow, multiple central kspace lines (beat beat) (KC)
- 20171114, collect inverse RECAR? at the same time as +RECAR  (KC)
- 20171201, detect abnormal flow during aquisition (somehow automatically track the center of kspace and flow position, then observe for beat difference... similar to navigator gating)?
