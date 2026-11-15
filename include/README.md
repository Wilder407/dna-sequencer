# Needleman-Wunsch Sequence Alignment in C++

A command-line tool that aligns two DNA or RNA sequences with the Needleman-Wunsch dynamic-programming algorithm. It reports the optimal alignment score, prints the alignment, and runs a shuffle test to show whether the score is higher than chance. Case studies: human **beta- vs. delta-globin** (HBB vs. HBD), a **sickle-cell point mutation** in HBB, and a reference **TFF1** mRNA vs. a TFF1 sequence from a breast-cancer cDNA library.

Lisa Wilder · originally built for CSCI 2270 (Data Structures), CU Boulder, spring 2025; refactored and extended in November 2026.

## Features

- Global alignment (classic Needleman-Wunsch) and an optional **semi-global** mode (`--free-end-gaps`) for aligning a partial read against a full sequence
- Configurable match, mismatch and gap scores
- Reads FASTA or plain-text files, or sequences typed or pasted at the prompt; strips whitespace, line numbers and headers
- Prints the alignment in 60-column blocks with a match line, plus matches, mismatches, internal and end gaps, and percent identity
- Prints the score and direction matrices for short sequences (useful for learning and checking by hand)
- **Shuffle test**: re-scores the alignment against N shuffled copies of one sequence, which keep its length and base composition, and reports how far the real score sits above them
- Unit tests: a textbook example, hand-checked cases, and 400 randomized checks that the traceback, the score and the memory-saving score function agree

## Quick start

```bash
make            # builds ./nw
make test       # builds and runs the unit tests
make demo       # runs all three case studies
```

Requires a C++17 compiler (g++ or clang++).

## Usage

```bash
./nw [options] [FILE1 [FILE2]]
```

| Input | Meaning |
|---|---|
| `FILE1 FILE2` | one sequence per file (FASTA or plain text) |
| `FILE1` | one file containing two records (headers start with `>` or `//`) |
| no files | type or paste each sequence, ending each with a blank line |

| Option | Default | Meaning |
|---|---|---|
| `--match N` | 1 | score added when two aligned bases are the same |
| `--mismatch N` | -1 | score added when two aligned bases differ |
| `--gap N` | -2 | penalty added for each base aligned to a gap |
| `--free-end-gaps` | off | gaps before the first or after the last base cost nothing |
| `--show-matrix` | auto | print the score and direction matrices (automatic when both sequences are 15 bases or fewer) |
| `--shuffles N` | 0 | compare the real score with N shuffled copies of the second sequence |
| `--seed N` | 42 | random seed for the shuffle test |

Example with the textbook sequences:

```text
$ printf 'GATTACA\n\nGCATGCU\n\n' | ./nw --gap -1
...
x  G-ATTACA
   | | |.|.
y  GCA-TGCU

Alignment score: 0
```

## How the algorithm works

Needleman-Wunsch finds the highest-scoring way to line up two whole sequences, allowing three kinds of columns:

- **Match or mismatch:** two bases are paired. Add `match` if they are the same, `mismatch` if not.
- **Gap in y:** a base of x is paired with a gap (`-`), meaning an insertion or deletion. Add the gap penalty.
- **Gap in x:** a base of y is paired with a gap. Add the gap penalty.

It fills a score matrix `F` of size (m+1) × (n+1), where `F[i][j]` is the best score for aligning the first `i` bases of x with the first `j` bases of y:

```
F[i][j] = max( F[i-1][j-1] + s(x_i, y_j),     ↖ pair x_i with y_j
               F[i-1][j]   + gap,             ↑ x_i against a gap
               F[i][j-1]   + gap )            ← y_j against a gap
```

Each cell depends on the three neighbors above and to the left:

| | |
|---|---|
| F[i-1][j-1] | F[i-1][j] |
| F[i][j-1] | **F[i][j]** |

The first row and column hold multiples of the gap penalty (0, −2, −4, …), because aligning `i` bases against nothing costs `i` gaps. A second matrix records which of the three moves gave the maximum.

**The alignment score is `F[m][n]`.** The traceback starts there and follows the stored directions back to `F[0][0]` to recover *which* bases pair up. It does not add anything to the score.

**What the gap penalty controls.** Gaps are only used when they pay for themselves. With match +1, mismatch −1 and gap −2, `ACGT` vs. `AGT` aligns as `ACGT / A-GT` for a score of 1 + (−2) + 1 + 1 = 1. A harsher penalty gives fewer gaps and more mismatches; a penalty near zero lets gaps appear everywhere.

**Global vs. semi-global.** Global alignment charges for gaps at the ends, which is right when both sequences cover the same region. When one sequence is a fragment or has extra material at an end (a primer, a poly-A tail, vector sequence), those end gaps are not meaningful, and `--free-end-gaps` sets the first row and column to 0 and takes the best score anywhere in the last row or column.

**Cost.** Filling the matrices takes O(mn) time and memory. The shuffle test only needs the score, so it uses a two-row version with O(n) memory. Aligning the 492- and 634-base TFF1 sequences with 200 shuffles takes about 0.3 seconds.

## Case study 1: beta-globin vs. delta-globin

HBB (beta-globin) and HBD (delta-globin) are neighboring human genes that arose from a gene duplication. Both coding sequences come from NCBI and are 444 bases long. As a check that the downloads were correct, the HBB sequence translates exactly to the beta-globin protein, and HBD differs from it at 10 amino acids, matching the known difference between the two proteins.

```bash
./nw data/HBB_NM_000518.5_cds.fasta data/HBD_NM_000519.4_cds.fasta --shuffles 1000
```

| Result (match +1, mismatch −1, gap −2) | |
|---|---|
| Alignment score | 381 |
| Identity | 93.0% (414 of 445 columns) |
| Mismatches / gap positions | 29 / 2 |
| Shuffled-score mean (SD), 1,000 shuffles | −60.5 (9.5) |
| Real score vs. shuffles | 46 SDs above the mean; 0 of 1,000 shuffles scored as high |

The genes are clearly related (93% identical, far above chance), and their 444-base sequences differ at 32 positions, which change 10 amino acids. The two gap positions sit next to each other and cancel out: with these scores, one pair of gaps costs 2 points less than the 3 mismatches it replaces. Both genes are the same length, so the biological picture is base substitutions, not insertions or deletions. Raising the gap penalty (`--gap -3`) removes the gaps and shows all 32 substitutions.

## Case study 2: detecting the sickle-cell mutation

Sickle-cell disease is caused by a single base change in HBB: c.20A>T, which turns codon 7 from GAG (glutamate) into GTG (valine). The file `data/HBB_sickle_c20A-T_cds.fasta` is the HBB coding sequence with that one base edited by hand.

```bash
./nw data/HBB_NM_000518.5_cds.fasta data/HBB_sickle_c20A-T_cds.fasta --shuffles 1000
```

```text
x  ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC
   |||||||||||||||||||.||||||||||||||||||||||||||||||||||||||||
y  ATGGTGCATCTGACTCCTGTGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAAC
```

The alignment scores 442 out of a possible 444, with **exactly one mismatch, at position 20**, and no gaps (99.8% identity). This is the case where pairwise alignment does find a point mutation: both sequences cover the same region, so any mismatch is a real difference.

## Case study 3: TFF1 in normal tissue vs. a breast-cancer library

**TFF1** (trefoil factor 1, also called pS2) is an estrogen-regulated gene that is often highly expressed in breast tumors. The question here: does the TFF1 sequence from a breast-cancer cDNA library differ from the reference?

| Sequence | Source | Length |
|---|---|---|
| x: reference TFF1 mRNA | NCBI [NM_003225](https://www.ncbi.nlm.nih.gov/nuccore/NM_003225) (as copied in 2025) | 492 bases |
| y: breast-cancer library TFF1 | NCBI [DN996463.1](https://www.ncbi.nlm.nih.gov/nuccore/DN996463.1) | 634 bases |

Both are in `data/m_RNASequences.txt`.

### Results (match +1, mismatch −1, gap −2, 200 shuffles)

| | Global | Free end gaps |
|---|---|---|
| Alignment score | 149 | 469 |
| Matches / mismatches | 479 / 2 | 475 / 6 |
| Gap positions (end / internal) | 6 / 158 | 164 / 0 |
| Shuffled-score mean (SD) | −159.4 (9.8) | 2.3 (1.6) |
| Shuffles scoring ≥ real score | 0 of 200 | 0 of 200 |

### What this shows

- **The sequences are the same transcript.** Bases 20–492 of the reference (473 bases, covering the full coding region and 3′ UTR) appear in the cancer-library sequence with no differences at all. The semi-global alignment has no internal gaps, and all 6 mismatches fall in the cancer-library sequence's first 8 bases.
- **All the differences are at the ends.** The cancer-library sequence starts at a different 5′ position and carries 153 extra bases at the 3′ end: a poly-A tail followed by sequence containing GCGGCCGC and TCTAGA, the NotI and XbaI restriction sites typical of cloning vectors. These are features of how the cDNA library was made, not changes in the gene.
- **Global alignment is the wrong tool for this pair.** It must account for the extra 3′ bases, so it scatters 158 gaps through the alignment and lowers the score to 149. Freeing the end gaps removes that artifact and raises the score to 469.
- **The similarity is far beyond chance** in both modes: no shuffled sequence came close to the real score.

So this comparison finds **no evidence of a mutation in TFF1**. That fits the biology: TFF1's link to breast cancer is mainly about how much of it is expressed, not changes to its sequence. Note that an early version of this README read the global score of 149 as "low similarity" pointing to a mutation. That was wrong: a raw score depends on sequence length and scoring, and only means something next to a baseline like the shuffle test.

## Limitations and next steps

- **Linear gap penalty.** Real insertions and deletions often span several bases; an affine penalty (a cost to open a gap plus a smaller cost to extend it) is more realistic.
- **Memory.** The full traceback needs O(mn) memory. Hirschberg's algorithm recovers the alignment in O(n) memory, which matters for sequences longer than a few thousand bases.
- **Bases only.** There is no support for ambiguity codes (N, R, Y…) or protein substitution matrices such as BLOSUM62.
- **One pairwise alignment is not mutation detection.** Detecting tumor mutations needs sequencing reads with quality scores, matched normal samples, and a variant-calling pipeline, not a single library sequence aligned to a reference.

## Project structure

```
include/Matrix.h                     data types and function declarations
src/Matrix.cpp                       matrix fill, traceback, input parsing, printing
app/main.cpp                         command-line program
tests/test_matrix.cpp                unit tests
data/HBB_NM_000518.5_cds.fasta       HBB coding sequence (NCBI)
data/HBD_NM_000519.4_cds.fasta       HBD coding sequence (NCBI)
data/HBB_sickle_c20A-T_cds.fasta     HBB with the sickle-cell c.20A>T change
data/m_RNASequences.txt              reference and breast-cancer-library TFF1 sequences
Makefile
```

## References

- Needleman, S. B. & Wunsch, C. D. (1970). A general method applicable to the search for similarities in the amino acid sequence of two proteins. *Journal of Molecular Biology*, 48(3), 443–453.
- Applied Bioinformatics, Chapter 3 (sequence alignment): https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-3/
- Interactive Needleman-Wunsch simulator, used to check results by hand: https://bioboot.github.io/bimm143_W20/class-material/nw/
- Biology 2e, *Cancer and Gene Regulation*, OpenStax (CC BY): https://courses.lumenlearning.com/suny-osbiology2e/chapter/cancer-and-gene-regulation/
- NCBI Nucleotide: [NM_000518.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5) (HBB), [NM_000519.4](https://www.ncbi.nlm.nih.gov/nuccore/NM_000519.4) (HBD), [NM_003225](https://www.ncbi.nlm.nih.gov/nuccore/NM_003225) (Homo sapiens TFF1 mRNA), [DN996463.1](https://www.ncbi.nlm.nih.gov/nuccore/DN996463.1) (breast-cancer library TFF1 sequence)

## Acknowledgments

The original version was written with help from ChatGPT-4o (OpenAI) in April–May 2025. The November 2026 refactor (semi-global mode, file input, shuffle test, tests, globin case studies and this README) was done with help from Claude (Anthropic).
