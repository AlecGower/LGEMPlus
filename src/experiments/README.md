# Experiments run for DS2023 paper submission

Here are the code files for recreating the experiments as submitted in the manuscript. There is code for the following experiments.

1. Initial and abduction of missing metabolites (not included in the ubiquitous or media) for:
    - yeastGEM;
    - iMM904; and
    - iND750.
2. Single gene deletion experiments.
    1. Using the logic model framework for each of the above three models.
    2. For yeastGEM with FBA evaluation.
    3. Evaluation of single gene deletion experiments with genome-wide screen data.
3. Abduction.
    1. Generation of hypotheses for false positive experiments from above.
    2. Extraction of pathways for each hypothesis.
    3. FBA constraint using pathways for each hypothesis.
    4. Single-gene deletion for each hypothesis.
    5. Ranking hypotheses based on above criteria.