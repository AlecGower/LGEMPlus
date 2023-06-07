# Experiments run for DS2023 paper submission

Here are the code files for recreating the experiments as submitted inthe manuscript. There is code for the following experiments.

1. Theory creation and single-gene deletion.
    1. Initial theory creation and abduction of missing metabolites (not included in the ubiquitous or media), and single-gene deletion experiments using the logic model framework for:
        - yeastGEM (Yeast8);
        - iMM904; and
        - iFF708.
    2. Single-gene deletion experiments using the logic model framework for the above models
    3. Single gene deletion experiments for yeastGEM with FBA evaluation.
    4. Evaluation of single gene deletion experiments with genome-wide screen data.
1. Abduction.
    1. Generation of hypotheses for false positive experiments from above.
    2. Extraction of pathways for each hypothesis.
    3. FBA constraint using pathways for each hypothesis.
    4. Single-gene deletion for each hypothesis.
    5. Ranking hypotheses based on above criteria.