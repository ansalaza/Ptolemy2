# Ptolemy2

Ptolemy2 is a stand-alone method for constructing gene-graphs given
long-read alignments to one or more reference genomes or pairwise whole-genome alignments.

The gene-graphs can then be used to analyze genome architecture--organization
and collection of genes in a genome--across a (metagenomic) sequencing dataset or a collection of genome assemblies.

This work is an extension of our previous published method, [Ptolemy](https://academic.oup.com/bioinformatics/article/34/17/i732/5093246),
and has thus far only been tested on genome assemblies and long-read metagenomic sequencing data of bacteriphages. But feel free to try it on on other organisms!

## How it works (in a nutshell...)

![](/figures/ptolemy2_overview.png?raw=true)

## Requirements
Java 1.8 or higher.

Optionally, [minimap2](https://github.com/lh3/minimap2) accessible from you PATH environment.

## How to run

Given long-read alignments in [PAF-format](https://github.com/lh3/miniasm/blob/master/PAF.md):
```bash
java -jar ptolemy2.jar build-db gene-graph -g all_gffs.txt \
--feature-types "CDS" \
--name-tag "product=" \
-r read_alignments.paf \
--prefix <some_string>
```

For whole-genome alignments, switch '-r' parameter with '-w'.

Output is a [GFA-formatted file](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md).

You can visualize the graph with [Bandage](http://rrwick.github.io/Bandage/) or use Ptolemy2's internal GFA-to-CSV converter and visualize with other
graph visualizers like [Gephi](https://gephi.org/).
