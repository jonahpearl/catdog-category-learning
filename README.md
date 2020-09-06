# catdog-category-learning
Doing the catdog analysis, the right way this time.

Note: the passive viewing code often refers to "trials" in the comments. I use "trial" as shorthand in my comments for "the presentation of an image". However, the Monkeys structure has a field named "TrialInfo" (which is itself a structure) which is not about individual image presentations. Instead, those are trials from the perspective of the behavioral task, which presents five images per trial. Most of this code doesn't deal with the TrialInfo structure -- it uses the CueInfo structure (which is intially constructed from the TrialInfo structure) to get all the times each image was presented and works from there.
(NB also, in the category training task, each trial really does only present one image.)
