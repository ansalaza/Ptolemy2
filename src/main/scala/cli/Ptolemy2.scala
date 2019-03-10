package cli

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
        "MAIN\n" +
        "graph-family          Construct gene-graph by inferring gene-families from pairwise gene alignments.\n"+
        "sequence-graph        Construct sequence-graph via MSA of defined protein families.\n" +
        "graph-alignment       Construct gene graph from long-read/whole-genome alignments.\n\n" +
        "SUPPLEMENT\n"+
        "path-metrics          Summarize metrics for all paths in a GFA file.\n" +
        "seq-db                Create a database for all sequences in given FASTA file.\n" +
        "gff2descriptions      Convert a list of GFF files into ptolemy2's descriptions file.\n" +
        "gfa-converter         Convert a GFA-file to various formats.\n\n" +
        "MISC. TOOLS\n" +
        "MINIMAP2 AUTOMATION TOOLS\n" +
        "index-genomes         Index a given list of genomes with minimap2\n" +
        "pairwise-alignment    Automate pariwsie genome alignment with minimap2\n\n" +
        "OTHERS\n" +
        "genbank-converter     Convert an assembly in genbank format to assembly (FASTA), proteins (FASTA), and " +
        "GFF3.\n" +
        "fetch-keys            Retain genes whose description matches a given list of string-keys\n" +
        "alignment-metrics     Obtain summary metrics for an alignment in PAF-format.\n" +
        "extract-genes         Extract gene-sequences from a genome assembly.\n"+
        "fetch-subsequences    Obtain (sub)sequences from a given read/assembly file\n\n"

      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "graph-family" => gene_graph.GeneGraphFamily.main(args.drop(1))
        case "sequence-graph" => gene_graph.SequenceGraph.main(args.drop(1))
        case "graph-alignment" => gene_graph.GeneGraphAlignment.main(args.drop(1))

        case "seq-db" => misc.Fasta2MapDB.main(args.drop(1))
        case "path-metrics" => metrics.ArchitectureMetrics.main(args.drop(1))
        case "gfa-converter" => misc.GFAconverter.main(args.drop(1))
        case "gff2descriptions" => misc.GFAconverter.main(args.drop(1))
        case "generate-labels" => misc.GenerateLabels.main(args.drop(1))

        case "index-genomes" => genome_alignment.IndexGenomes.main(args.drop(1))
        case "pairwise-alignment" => genome_alignment.PairwiseAlignment.main(args.drop(1))

        case "genbank-converter" => misc.GenbankConverter.main(args.drop(1))
        case "fetch-keys" => misc.FetchWithStrings.main(args.drop(1))
        case "alignment-metrics" => metrics.AlignmentMetrics.main(args.drop(1))
        case "extract-genes" => misc.ExtractGenes.main(args.drop(1))
        case "fetch-subsequences" => misc.FetchSubSequences.main(args.drop(1))

        case _ => println(help)
      }
    }
  }

}
