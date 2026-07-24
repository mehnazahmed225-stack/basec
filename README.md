## ExonScorer

ExonScorer is a transcript-to-genome alignment algorithm that estimates the level of transcript support for each nucleotide in a genomic sequence. The underlying assumption is that mature mRNA molecules are derived from expressed exons, while intronic regions are removed during RNA splicing.

Unlike traditional exact-match approaches, ExonScorer searches for **ordered subsequence matches** between an RNA transcript and genomic DNA. This allows transcript sequences to align across intervening genomic regions while preserving nucleotide order.

### Example

RNA sequence:

```text
ATC
```

Genomic sequence:

```text
AGTGCG
^ ^ ^
```

Although `ATC` does not appear as a contiguous substring, it exists as an ordered subsequence of the genomic sequence and is therefore considered a valid transcript-supported path.

### Algorithm

For each genomic suffix:

1. The suffix is considered only if it begins with the same nucleotide as the RNA sequence.
2. The algorithm determines whether the RNA sequence can be embedded within the suffix as an ordered subsequence.
3. If a valid subsequence exists, the positions participating in the match are recorded.
4. A compactness score is calculated:

$$
Weight = \frac{RNA\ Length}{Match\ Span}
$$

where:

```text
Match Span = Last Matched Position - First Matched Position + 1
```

This scoring scheme rewards compact alignments while penalizing highly fragmented alignments.

### Examples

Perfectly compact match:

```text
RNA: ATCG
DNA: ATCG

Weight = 4 / 4 = 1.0
```

Fragmented match:

```text
RNA: ATCG
DNA: A.....T.....C.....G

Weight = 4 / 16 = 0.25
```

As the number of skipped nucleotides increases, the contribution of the match decreases.

### Exon Proclivity Calculation

Each nucleotide participating in a valid transcript alignment receives the calculated weight as evidence of exon membership.

For every genomic position, all supporting weights are accumulated and averaged:

```text
Exon Proclivity =
Average(Transcript Support Weights)
```

Positions that repeatedly participate in compact transcript alignments receive higher scores, while positions that rarely participate or only participate in highly fragmented alignments receive lower scores.

### Interpretation

ExonScorer generates a per-nucleotide **Exon Proclivity Score**, which can be visualized as a genomic support track or heatmap.

Higher scores indicate:

- Strong transcript support
- Frequent participation in valid RNA-DNA alignments
- Compact exon-like alignment patterns

Lower scores indicate:

- Weak transcript support
- Sparse participation in transcript alignments
- Highly fragmented or absent matches

### Current Limitations

- Only perfect nucleotide matches are currently supported.
- Insertions, deletions, and mismatches are not yet modeled.
- Only the first valid subsequence path is evaluated for each suffix.
- The Exon Proclivity Score should be interpreted as a transcript-support metric rather than a formal probability of exon membership.

### Future Development

Planned improvements include:

- Mismatch-tolerant alignments
- Gap penalties based on intron-length distributions
- Multiple-path subsequence scoring
- Splice-site signal integration
- Machine-learning-based exon probability prediction

The long-term goal of ExonScorer is to provide a lightweight framework for identifying transcript-supported genomic regions and generating exon evidence profiles that can be integrated into downstream gene annotation pipelines.
