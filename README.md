fMRI data + stimuli loading for SOC project
==========================================

This respository contains the code used to generate the stimuli
presented in the Winawerlab Second-Order Contrast project experiments
from 2015-2016, as well as code to load the data from those
experiments in a canonical format with plenty of metadata.

First, make sure to point the `soc_fmripath.m` and `soc_stimpath.m`
scripts at the appropriate locations to find the fMRI data and the
stimuli. These are by default set to point to the winawerlab server,
but you can point them at a local copy of the data if you want.


Experiments
===========

This repository pertains to four individual fMRI subject datasets:

  * 2015-04-07 - subj004

    * This pilot study was run just before VSS, to test two new image
      types - "waves" and "noise bars" - alongside two existing image
      types - "gratings" and "patterns". Only these four stimulus
      types were run.

  * 2015-06-19 - subj022

    * This is the first "full" dataset collected at NYU for the SOC
      experiment. A variety of images were used: patterns, gratings,
      noisebars, and waves, at different sparsities, orientations, and
      contrasts; as well as "cross-modulated" gratings. No retinotopic
      stimuli were used, becasue the hope was that a separate
      retinotopy scan for this subject could be used. This dataset is
      designated the "horizontal" dataset, because the default
      orientation for gratings was horizontal.

  * 2015-10-22 - subj001

    * This is the second "full" dataset collected at NYU for the SOC
      experiment. The stimuli used are basically identical to those
      used in June, except this is the "vertical" dataset, in that
      every stimulus was rotated 90 degrees before display, making the
      default orientation consist of vertical bars. The
      gratings-versus-patterns discrepancy is less pronounced at
      vertical orientations, suggesting some orientation-dependence.

  * 2016-02-17 - subj001

    * This small dataset explores spatial frequency tuning, with the
      intent to inform the design of a multiband model. Pattern
      stimuli at different spatial frequencies were used. Noise
      stimuli were generated, but not actually used in the experiment.


Code organization
=================
The subdirectories of this repository are:
  
  * dataload (`load_subj` scripts)
  
    * These scripts load in all the data and stimuli from a given run
      into a "data" struct which contains metadata
      
  * imagereshape
  
    * Stacks of images are represented in different ways across this
      code and KNK's code. Stimuli consist of images (X * Y pixels) of
      different categories (C) where each category can have multiple
      examples (frames, F), and each individual image can be
      represented in a multi-band format with several bands (B).

        * cell: A cell array of C categories, each cell containing a
          matrix of images shaped X pixels * Y pixels * F frames. Used
          in the file `stimulus.mat` downloadable from KNK's
          website, but nowhere else in the code.
				      
	* stack: A matrix X * Y * C * F. Can also have a fifth
          dimension B for bands. This is the structure that most
          easily allows slicing and indexing, for visualizing images.
  
        * flat: A matrix (X * Y) * (C * F) * B bands, where each image has
	  been turned into a vector, and category boundaries have been
	  removed, although band distinctions remain. The "model
	  steps" mostly operate on flattened images for efficiency
	  purposes
										    
	* pixel vector: A matrix (X * Y) * C * F, where each image has
          been turned into a vector, but category boundaries
          remain. KNK's model-fitting code expects something like a
          pixel vector, but with the dimensions flipped to be C *
          (X * Y) * F, so interfaces with KNK's code will often include
          a "permute" step

  * knkcode
  
     * All the code in this repository was originally written with a
       dependency on the `knkutils` package. For explicitness, the
       specific dependencies have been pulled out.


  * stimulusgen

     * `stimulusgen_<date>` scripts contain exactly what was run to
        generate the stimuli shown for a specific experiment.
        Although the random seed was unfortunately not always
        captured, the files saved as `stimuli-<date>.mat` /should/
        contain the actual original stimuli used, whereas
        `stimuli-<date>_rerun.mat` were regenerated at a later point
        to test the script. However, if in doubt, you can always
	go fetch the .mat file from the "Stimulus" folder of the
	fMRI_CBI folder, which is guaranteed to be the one that was
	actually shown to the subject.
	
     * The other functions in this directory are helper functions to
       generate certain classes of stimuli, bandpass filters, and so
       on.

  * utils

     * These functions capture a few common motifs in masking and indexing.
       
										    


