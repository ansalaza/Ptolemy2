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
        "gene-graph          Construct a gene graph from a given set of genomes.\n\n" +
        "MINIMAP2 AUTOMATION TOOLS\n" +
        "index-genomes         Index a given list of genomes with minimap2\n" +
        "pairwise-alignment    Automate pariwsie genome alignment with minimap2\n\n" +
        "MISC. GRAPH TOOLS\n" +
        "gfa-converter         Convert a GFA-file to various formats.\n"
      )

    if (args.length == 0) {
      println(help)
    } else {
      args(0) match {
        case "gene-graph" => gene_graph.GeneGraph.main(args.drop(1))

        case "index-genomes" => genome_alignment.IndexGenomes.main(args.drop(1))
        case "pairwise-alignment" => genome_alignment.PairwiseAlignment.main(args.drop(1))

        case "gfa-converter" => gene_graph.GFAconverter.main(args.drop(1))
        case _ => println(help)
      }
    }
  }

}
