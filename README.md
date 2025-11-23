# RunningExercise2_BINP16
### Oliver Todreas, Oct 21 2025.
This project contains programs, inputs, outputs, tests, and documentation. This includes:

```markdown
FastaAlignerProject-1/
├── FastaAligner.py
├── FastaAlignerPlotter.py
├── README.md
└── testdata/
    ├── score.fna
    ├── parameters.txt
    ├── failtests/
    |   ├── incorrect_fasta_type.txt
    |   ├── corrupted_fasta_file.fna
    |   ├── invalid_characters_warning.fna
    |   ├── single_sequence_error_fasta.fna
    |   ├── unequal_seq_lengths_fasta.fna
    |   ├── incorrect_param_type.fna
    |   ├── misconfigured_param_noequals.txt
    |   ├── misconfigured_param_multipleequals.txt
    |   ├── misconfigured_param_paramnames.txt
    |   └── misconfigured_param_nonnumerical.txt
    └── outputs/
        ├── outputs_scorefna.txt
        └── dotplots_scorefna/
            └── id1_id2.png
```

Extensive documentation is provided inside the python files. Note that user-defined functions have
their own documentation.

To test for errors, run the programs as prompted to in the program documentation, but pass
testdata/failtests/ files as inputs.

Passing testdata/score.fna as the input to both programs will return identical outputs to those
provided in testdata/outputs/. To get identical results, do not pass a parameters.txt file to
FastaAligner.py.