### Contents of *ketaconf/expe* subfolder

This subfolder contains the core functions and scripts used run the KETACONF experiment.

**KETACONF_runner.m** -- Runner script for launching the KETACONF experiment. <br />
This is a runner script for launching the KETACONF experiment.

**KETACONF_run_expe.m** -- Run KETACONF experiment. <br />
This is the main function for running the KETACONF experiment.

**KETACONF_gen_expe.m** -- Generate KETACONF experiment. <br />
This function is called by the main function that runs the experiment.

**gen_seqang.m** -- Generate sequence angles with titrated evidence.

**make_card.m** -- Generate card stimulus.

**make_pie.m** -- Generate pie chart stimulus for lottery.

**update_sd_cue.m** -- Update effective inference noise s.d. from provided responses. <br />
This function is used by the online staircasing procedure to adjust task difficulty throughout the experiment.

**itrnd.m** -- Random number generator using inverse transform method.

**vmk2r.m** -- Get von Mises coherence from concentration.

**vmr2k.m** -- Get von Mises concentration from coherence.

**vmpar.m** -- Get von Mises parameters from array of angles.

The *expe/Documents* subfolder contains the instructions for the KETACONF experiment, written in French.

The *expe/Toolboxes* subfolder contains additional low-level functions that are used to run the KETACONF experiment (input-output functions that build on core Psychtoolbox-3 functions, random-number generation functions, psychophysical staircasing functions, and stimulus generation and presentation functions).
