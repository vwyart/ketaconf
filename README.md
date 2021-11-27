# KETACONF
Core functions and scripts for KETACONF study (Salvador et al., 2022, *Nature Comm.*)

Because tested participants did not provide explicit written consent regarding the posting of their anonymized data on public repositories, individual anonymized behavioral and EEG datasets are available upon request at *valentin.wyart@inserm.fr*.

The experiment code for running the study is available in the *ketaconf/expe* subfolder. The code runs using MATLAB and the Psychtoolbox-3 toolbox for MATLAB (*http://psychtoolbox.org/*). The analysis code for fitting the inference and confidence models to behavior, for analyzing the EEG data, and for simulating the premature commitment model of NMDA receptor hypofunction, is available in the *ketaconf/analysis* subfolder.

### Contents of *ketaconf/analysis* subfolder

**classify_proj.m** -- Decode discrete labels using simple linear classifier based on multivariate normal model.

**decode_pinv.m** -- Decode continuous variable from multivariate data using pseudoinverse (i.e., an inverted encoding model).

**decode_pinv_tgen.m** -- Decode continuous variable from multivariate data using pseudoinverse (i.e., an inverted encoding model) with temporal generalization.

**fit_model_resp.m** -- Fit noisy inference model to binary responses.

**fit_model_noisesource.m** -- Fit inference model with different noise sources.

**fit_model_conf.m** -- Fit noisy inference and confidence model to binary responses and opt-out judgments.

**get_aroc.m** -- Get area under ROC curve.

**get_proc.m** -- Get multinomial ROC sensitivity to categorical labels.

**get_sesstype.m** -- Get session type for KETACONF study.

**load_subjsess.m** -- Load subject/session behavioral data for KETACONF study. <br />
*This script requires individual anonymized behavioral datafiles to run.*

**mvnllh.m** -- Multivariate normal log-likelihood.

**normllh.m** -- Normal log-likelihood.

**read_dat_mod.m** -- Read stimulus-locked EEG data of KETACONF study. <br />
*This script requires individual anonymized behavioral and EEG datafiles to run.*

**read_dat_resp.m** -- Read response-locked EEG data of KETACONF study. <br />
*This script requires individual anonymized behavioral and EEG datafiles to run.*

**script_poweranalysis.m** -- Analysis script for response preparation power analysis. <br />
*This script requires individual anonymized behavioral and EEG datafiles to run.*

**script_simmodel_precom.m** -- Run simulations of the pre-commitment model of inference and confidence.

**script_symptomscores.m** -- Analyze symptom scores and ketamine plasma concentration for KETACONF study. <br />
*This script requires individual the anonymized symptom scores datafile to run.*

### Contents of *ketaconf/expe* subfolder

**KETACONF_runner.m** Runner script for launching the KETACONF experiment. <br />
This is a runner script for launching the KETACONF experiment.

**KETACONF_run_expe.m** Run KETACONF experiment. <br />
This is the main function for running the KETACONF experiment.

**KETACONF_gen_expe.m** Generate KETACONF experiment. <br />
This function is called by the main function that runs the experiment.

The *expe/Documents* subfolder contains the instructions for the KETACONF experiment, written in French.

The *expe/Toolboxes* subfolder contains additional low-level functions that are used to run the KETACONF experiment (input-output functions that build on core Psychtoolbox-3 functions, random-number generation functions, psychophysical staircasing functions, and stimulus generation and presentation functions).
