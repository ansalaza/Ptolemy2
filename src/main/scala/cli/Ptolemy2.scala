package cli

import gene_graph.GeneGraphFamily
import misc.{ExtractGenes, GFAconverter}

/**
  * Author: Alex N. Salazar
  * Created on 7-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object Ptolemy2 {

  def main(args: Array[String]): Unit = {
    val help = (
      "Usage: java -jar ptolemy2.jar [tool]\n\n" +
        "PTOLEMY2 (MAIN)\n" +
        "graph-family          Construct gene-graph by inferring gene-families from pairwise gene alignments.\n"+
        "graph-alignment       Construct gene graph from long-read/whole-genome alignments.\n" +
        "path-metrics          Summarize metrics for all paths in a GFA file.\n" +
        "gfa-converter         Convert a GFA-file to various formats.\n\n" +
        "PTOLEMy2 (SUPPLEMENT)\n"+
        "gff2descriptions      Convert a list of GFF files into ptolemy2's descriptions file.\n\n" +
        "MISC. TOOLS\n" +
        "MINIMAP2 AUTOMATION TOOLS\n" +
        "index-genomes         Index a given list of genomes with minimap2\n" +
        "pairwise-alignment    Automate pariwsie genome alignment with minimap2\n\n" +
        "OTHERS\n" +
        "genbank-converter     Convert an assembly in genbank format to assembly (FASTA), proteins (FASTA), and " +
        "GFF3.\n" +
        "fetch-keys            Retain genes whose description matches a given list of string-keys\n" +
        "alignment-metrics     Obtain summary metrics for an alignment in PAF-format.\n" +
        "extract-genes         Extract genes.\n"+
        "fetch-subsequences    Obtain (sub)sequences from a given read/assembly file\n\n"

      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "graph-family" => GeneGraphFamily.main(args.drop(1))
        case "graph-alignment" => gene_graph.GeneGraphAlignment.main(args.drop(1))
        case "path-metrics" => metrics.ArchitectureMetrics.main(args.drop(1))
        case "gfa-converter" => GFAconverter.main(args.drop(1))

        case "gff2descriptions" => GFAconverter.main(args.drop(1))

        case "index-genomes" => genome_alignment.IndexGenomes.main(args.drop(1))
        case "pairwise-alignment" => genome_alignment.PairwiseAlignment.main(args.drop(1))

        case "genbank-converter" => misc.GenbankConverter.main(args.drop(1))
        case "fetch-keys" => misc.FetchWithStrings.main(args.drop(1))
        case "alignment-metrics" => metrics.AlignmentMetrics.main(args.drop(1))
        case "extract-genes" => ExtractGenes.main(args.drop(1))
        case "fetch-subsequences" => misc.FetchSubSequences.main(args.drop(1))

        case _ => println(help)
      }
    }
  }

}
