### Contents of *ketaconf/analysis* subfolder

This subfolder contains the core functions and scripts used to analyze the behavior and EEG data of the KETACONF study.

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
