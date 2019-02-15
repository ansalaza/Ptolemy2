package utilities

import utilities.GFAutils.{GFAreader, GFAwriter}
import utilities.GFFutils.Gene
import utilities.GeneProjectionUtils.FingerTree

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 14-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GeneGraphUtils extends GFAwriter with GFAreader {

  /**
    * Type alias for gene graph
    */
  type GeneGraph = Map[Int, List[Int]]

  /**
    * Empty gene graph
    */
  val empty_gene_graph = Map.empty[Int, List[Int]]

  /**
    * Type alias for paths
    */
  type Paths = Map[String, List[Gene]]

  /**
    * Empty paths
    */
  val empty_paths = Map.empty[String, List[Gene]]

  /**
    * Type alias for node coverage
    */
  type NodeCoverage = Map[Int, Int]

  /**
    * Empty node coverage
    */
  val empty_node_coverage = Map.empty[Int, Int]

  /**
    * Type alias for edge coverage
    */
  type EdgeCoverage = Map[(Int, Int), Int]

  /**
    * Empty edge coverage
    */
  val empty_edge_coverage = Map.empty[(Int, Int), Int]

  /**
    * Method to identify all connected components in a given graph
    *
    * @param graph Graph as a Map(node ID -> List[Node ID]) (Map[Int, List[Int])
    * @return Connected components as list of sets of node IDs (List[Set[Int])
    */
  def getConnectedComponents(graph: Map[Int, List[Int]]): List[Set[Int]] = {

    /**
      * Tail-recusrive method to traverse through a given list of nodes and identify the the corresponding connected
      * component
      *
      * @param toVisit List of node IDs to traverse
      * @param visited Set of node IDs already visited
      * @return Set of node IDs visited during the traversal (e.g. connected component)
      */
    @tailrec def breadthFirstTraversal(toVisit: List[Int], visited: Set[Int]): Set[Int] = {
      toVisit match {
        //no more remaining nodes to visit
        case Nil => visited
        //get next node to visit
        case node :: tail => {
          //skip if current node is already visited
          if (visited(node)) breadthFirstTraversal(tail, visited)
          //add neighbouring nodes to node to visit and and current node to set of nodes already visited
          else breadthFirstTraversal(graph.getOrElse(node, List[Int]()) ::: toVisit, visited + (node))
        }
      }
    }

    /**
      * Tail-recusrive method to perform a breadth first traversal through all nodes in the original graph provided
      * and identify all connected components
      *
      * @param remaining_nodes Remaining nodes to traverse
      * @param ccs             Accumulating connected components
      * @return All connected components
      */
    @tailrec def _getConnectedComponents(remaining_nodes: List[Int],
                                         ccs: List[Set[Int]]): List[Set[Int]] = {
      //node more nodes to visit
      remaining_nodes match {
        //no more nodes to visit
        case Nil => ccs
        //traverse head node
        case head :: tail => {
          //perform breadth-first traversal starting from the given head node to get the local connected component
          val local_cc = breadthFirstTraversal(List(remaining_nodes.head), Set())
          //remove all nodes in local cc in the remainindg nodes and move on
          _getConnectedComponents(remaining_nodes.filterNot(local_cc(_)), local_cc :: ccs)
        }
      }
    }

    //get all nodes in the given graph
    val all_nodes = graph.foldLeft(List[Int]())((list, entry) => list ::: (entry._1 :: entry._2))
    //identify all connected components
    _getConnectedComponents(all_nodes, List())
  }


  /**
    * Method to construct native gene graph given a fingertree data structure
    *
    * @param fingertree    Finger tree of some genome
    * @param initial_graph optional gene graph to build-on from
    * @return GeneGraph
    */
  def fingertree2GeneGraph(fingertree: FingerTree,
                           initial_graph: GeneGraph = empty_gene_graph): (GeneGraph, List[Gene]) = {
    //get path by gene ids and sort by id
    val path = fingertree.toList.sortBy(_._1._1).map(_._2)
    //iterate through each edge and create gene graph
    val updated_graph = path.sliding(2).foldLeft(initial_graph)((acc_graph, genes) => {
      //get nodes
      val (node1, node2) = (genes.head.id, genes(1).id)
      //get edges
      val node1_edges = node2 :: acc_graph.getOrElse(node1, List[Int]())
      val node2_edges = node1 :: acc_graph.getOrElse(node2, List[Int]())
      //update graph
      (acc_graph + (node1 -> node1_edges)) + (node2 -> node2_edges)
    })
    (updated_graph, path)
  }

  /**
    * Function to create a non-directional syntenic-anchor graph based on a list of 2-tuples representing syntenic
    * anchors. Returns map where key is gene ID and values is a list of IDs representing neighbouring gene IDs
    * @return Map[Int, List[Int]
    */
  def createSAgraph: List[(Gene, Gene)] => Map[Int, List[Int]] = syntenic_anchors => {
    //join map of gene -> edges and edges -> gene
    (syntenic_anchors.groupBy(_._1).mapValues(_.map(_._2)) ++ syntenic_anchors.groupBy(_._2).mapValues(_.map(_._1)))
      //create map of INTs
      .map(x => (x._1.id, x._2.map(_.id)))
  }

  /**
    * Define canonical (new) IDs given a list of connected components as list of set of gene IDs.
    * @return Map[Int,Int]
    */
  def defineCanonicalIDs: List[Set[Int]] => Map[Int,Int] = sa_ccs => {
    //iterate through each connected component
    sa_ccs.foldLeft(List[(Int, Int)]())((gene_alignments, cc) => {
      //get lowest gene ID
      val merged_gene_id = cc.min
      //iterate through each gene and add tuple making old gene id -> new gene id
      cc.toList.foldLeft(gene_alignments)((local_gene_alignments, gene) => {
        (gene, merged_gene_id) :: local_gene_alignments
      })
    }).toMap
  }
}
