## GRITIC Scripts Repository
This is the repository for the analysis scripts accompanying the manuscript for *The history of chromosomal instability in genome doubled tumors*. The repository for GRITIC, the tool developed in this work can be found [here](https://github.com/VanLoo-lab/gritic).

## GRITIC  scripts
These scripts are used run GRITIC on the PCAWG and Hartwig cohorts. They are also used to build simulated samples used to evaluate GRITIC, using PCAWG and Hartwig samples as the initial framework. Access to the raw PCAWG and Hartwig data is required to run these scripts.

- ```HartwigDataLoader.py``` - A data loader for the Hartwig cohort.
- ```hartwig_handler.py``` Runs GRITIC on the Hartwig cohort.
- ```PCAWGDataLoader.py``` - A data loader for the PCAWG cohort..
- ```hartwig_handler.py``` Runs GRITIC on the PCAWG cohort.
- ```SampleSimulator.py``` Produces simulated tumor samples from a template.
- ```state_validation_run.py``` Simulates samples and runs through GRITIC.
- ```state_validation_run_titrate.py``` Simulates samples and runs through GRITIC.
## Analysis scripts
These scripts produce the figures for the manuscript.
- ```arm_pre_post.py``` - A precursor script to generate a table for figures 4A,B and supplementary figures 37-41. Requires ```pre_post_non_gain.py``` and ```pre_post_non_loss.py``` to be run first.
- ```cancer_type_parsimony.py``` Produces figure 2D and supplementary figures 15D,E.
- ```cn_wgd_figures.py``` Produces figures 1A-C and supplementary figure 1.
- ```combinesamplegainplots.py``` Produces supplementary figures 19 and 26.
- ```DataTools.py``` Helper functions for data loading.
- ```evaluate_parsimony_penalty.py``` Evalutes the effect of the non-parsimony penalty on parsimony inference, a precursor to supplementary figures 17 and 18.
- ```event_proportions_relative_to_wgd_probabilistic.py``` Produces figures 4C,D and supplementary figures 42-49.
- ```fraction_of_gains_job.py``` Produces figure 3D and supplementary figures 28A-B. Requires ```pre_post_permute_job.py``` to be run first.
- ```gain_distributions_across_samples.py``` Produces figure 3C.
- ```GainClassificationAnalysis.R``` Produces supplementary figure 36, requires access to raw PCAWG data.
- ```GainClassifier.R``` Helper script for ```GainClassificationAnalysis.R```.
- ```generic_gain_histograms.py``` Produces supplementary figure 27.
- ```hartwig_multiregion.py``` Produces supplementary figure 12. Requires raw data to run.
- ```hartwig_pcawg_copy_number_export.R``` A precursor script for ```cn_wgd_figures.py```.
- ```major_3_4_analysis.py``` Produces supplementary figure 20.
- ```major_cn_timing_job.py``` Produces figures 2E,F and supplementary figure 22.
- ```make_sample_timing_plots.py``` Plots sample gain timing posteriors, a precursor to supplementary figures 19 and 26.
- ```pan_gain_probabilistic.py``` Produces figure 3B and supplementary figures 24 and 25.
- ```parsimony_plots_with_bootstrapping.py``` Produces figure 3D and supplementary figures 15A,B and 16.
- ```plot_parsimony_prior_evaluations.py``` Produces supplementary figures 17 and 18.
- ```plot_wgd_sampling.py``` A helper script to produce figure 3B and supplementary figures 24 and 25.
- ```plot_wgd_timing_event_proportions.py``` Produces figures 3F,G and supplementary figures 30-33.
- ```plot_passage_cn_changes.R``` Produces figure 2B. Requires ```process_cn_passage_data.py``` to be run first.
- ```pre_post_correlations_corrected.R``` Produces figures 4A,B and supplementary figures 37-41. Requires ```arm_pre_post.py``` to be run first.
- ```pre_post_major_cn.py``` Produces supplementary figures 21.
- ```pre_post_non_gain.py``` Precursor script to produce aggregated gain timing files. Required for ```arm_pre_post.py``` and ```event_proportions_relative_to_wgd_probabilistic.py```.
- ```pre_post_non_loss.py``` Precursor script to produce aggregated loss timing files. Required for ```arm_pre_post.py``` and ```event_proportions_relative_to_wgd_probabilistic.py```.
- ```pre_post_permute_job.py``` Precursor script to produce aggregated loss timing files.
- ```process_cn_passage_data.py``` Precursor script to produce figure 2B
- ```route_calibration_plot.py``` Produces supplementary figures 8-11.
- ```route_difference_analysis.py```Produces supplementary figures 13-14.
- ```timing_matchup_sim.py``` Produces figures 1H and supplementary figures 3-7.
-  ```usarc_analysis_medicc.ipynb``` Produces figure 3A.
- ```wgd_calling_analysis_script.py``` Produces supplementary figure 51.
- ```wgd_constraint_analysis.py``` Produces supplementary figure 2.
- ```wgd_dist_permutation_job.py``` Produces figure 3E and supplementary figures 28C,D and 29.
