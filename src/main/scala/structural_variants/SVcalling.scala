package structural_variants

import java.io.File
import utilities.SequenceGraphUtils.msa2SequenceGraph
import utilities.GFAutils.GeneGraphReader

/**
  * Author: Alex N. Salazar
  * Created on 10-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object SVcalling extends GeneGraphReader {

  def main(args: Array[String]): Unit = {
    val gfa = new File("D:\\asalazar\\Desktop\\playground\\phage_architecture\\t1\\gene_graph.gfa")
    val (graph,paths) = loadGFA(gfa)

    paths2EdgeLabels(paths).filter(x => x._2("EC159") && x._2("EC131_SA29")).toList.sortBy(_._2.size).foreach(println)

    val msa = new File("D:\\asalazar\\Desktop\\playground\\phage_architecture\\t1\\test.mfsa.fa")

    val (seqgraph, seqpaths) = msa2SequenceGraph(msa)

    //seqgraph.toList.sortBy(_._1.id).foreach(println)
    //seqpaths.foreach(println)


  }

}
