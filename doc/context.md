# Project Context: Viral Mimicry & Host-Range Assignment via Structural Search

> Paste this into a new chat to give Claude the background of the ongoing project.
> It captures established understanding, open decisions, and important caveats.

## What I'm doing
I'm using structural search (Foldseek) to detect **viral protein structural mimicry of
host proteins**, with the goal of using that mimicry signal to say something about a
virus's **host range**. I query viral proteomes against host proteomes / a labeled target
DB and look for structural neighbors.

### Settled decisions (current phase)
- **I only care about STRUCTURAL mimicry, not functional mimicry.** Fold match is the
  signal I want; whether it's functionally real is out of scope for now.
- **First goal = reproduce Lasso at the broad TAXONOMIC-DIVISION level** (bacteria /
  plant-fungi / invertebrate / vertebrate), then move to higher resolution later.
- **I'm fine using pre-annotated hosts** (virus-hostDB style), like Lasso did. Not trying
  to predict host range from scratch yet — reproducing the correlation first.
- **I will handle proteome-size/composition bias with a statistical CORRECTION
  (hypergeometric), NOT by matching proteome sizes.** No size-matching hassle.
- **Query side:** virus taxid -> BFVD -> download that virus's predicted protein
  structures -> use as Foldseek query input. (BFVD = predicted viral structures,
  ColabFold/AF2-based.)

## What I want or considering

### Target database decision (division-level reproduction)

**As a starting point I am considering using Foldseek's prebuilt `Alphafold/Proteome` as the target DB**
- It's the AFDB "Proteome" subset: ~564,000 structures = the human proteome + other
  key organisms. Spans all my divisions:
  bacteria, fungi, plants, invertebrates (incl. Aedes aegypti mosquito, which Lasso had
  to model by hand), vertebrates incl. human, plus some protozoan parasites.
- Small: fits in <4 GB RAM, fast to search

** Feel free to let me know if you think there is a better idea for what target DB to use (e.g. the entire UniProt reference proteomes set ~36,465 proteomes, AlphaFold/UniProt50, etc)**

**Later (higher resolution):** swap `Alphafold/Proteome` or whatever we end up using for a custom set of per-genus or per-level reference proteomes (pull from AFDB by proteome ID, or a taxonomy-labeled slice of AFDB/UniProt50) and bin at a finer taxonomic level. But again, feel free to suggest other better ideas.

### Similarity cutoff
Lasso used a conservative GLOBAL criterion (SAS < 2.5 Angstrom) to avoid counting local/
partial similarity as mimicry. Foldseek defaults are more local + E-value driven, so out
of the box I'll get a looser, larger hit set than he did.
- To approximate his global-fold definition: I could run Foldseek in **TM-align mode
  (`--alignment-type 1`)** and threshold on **TM-score (~0.5 = "same fold")**, and/or
  require high bidirectional coverage (`-c 0.7-0.8`, `--cov-mode 0`) with a strict E-value.
- Without this, mimicry counts are inflated vs Lasso and enrichment may look different even
  if the pipeline is otherwise correct.

## What I want you to do

1. Using what you know about Lasso's paper and the context above, propose a concrete plan for you would reproduce Lasso's results using BFVD and Foldseek. Give me a step-by-step pipeline plan.

2. Suggest how you would choose the target database (e.g. the entire UniProt reference proteomes set ~36,465 proteomes, AlphaFold/UniProt50, etc), or if you aim to create a custom one, how you would do it.

3. Suggest how you would choose the similarity cutoff (e.g. TM-score, coverage, E-value, etc), or if you aim to create a custom formula, how you would do it.


