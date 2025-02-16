# catdog-category-learning
Code for working with pre-processed data from the paper [Neural Correlates of Category Learning in Monkey Inferior Temporal Cortex](https://doi.org/10.1523/JNEUROSCI.0312-24.2024). 

This code was written + run with Matlab 2020a, though should be largely forward-compatible. Additional dependencies can be downloaded from: https://github.com/jonahpearl/matlab_utils. Pre-processed data is available for download from Zenodo at: https://zenodo.org/records/14787249.

To reproduce figures, first download the data from the Zenodo page above. Then, go through each of the files labeled `figs_figX...` and either run the script or follow the pointers to other scripts. Please file an issue here with any questions.

Note: the passive viewing code (scripts with the "pv" prefix) often refers to "trials" in the comments. I use "trial" as shorthand in my comments for "the presentation of an image". However, the Monkeys structure contains a structure named "TrialInfo" which is not about individual image presentations. Instead, those are trials from the perspective of the behavioral task, which presents five images per trial. Because most of this code doesn't deal with the TrialInfo structure -- it uses the CueInfo structure (which is intially constructed from the TrialInfo structure) to get all the times each image was presented and works from there -- this shouldn't be a big problem. Just FYI in case someone looks closely.
(NB also, in the category training task (scripts with the "skip" prefix), each trial really does only present one image, on which the monkey must release the bar in one of two intervals to classify the image as a cat or dog.)
