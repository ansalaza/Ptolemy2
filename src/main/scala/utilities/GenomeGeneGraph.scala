package utilities

import utilities.AlignmentUtils.Alignment
import utilities.FileHandling.timeStamp
import utilities.GFFutils.Gene
import utilities.GeneProjectionUtils.{FingerTree, filterByCoverage, normalizePosition}
import utilities.GeneGraphUtils._

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 21-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GenomeGeneGraph {

  /**
    * Tail-recursive method to find syntenic anchors based on whole genome alignment. Pseudocode as follows:
    *
    * Initialize empty list of 2-tuples representing 'syntenic anchors'
    * For each genome:
    * For each alignment of current genome in sorted order:
    * Get all overlapping genes in query genome based on current alignment coordinates
    * For each overlapping gene in the query:
    * Project coordinates of overlapping gene onto target genome
    * Identify overlapping gene(s) in target genome based on projected coordinates
    * For each overlapping target gene:
    * Add (query gene, target gene) to list of syntenic anchors
    *
    * @param alignments Curated list of all alignments per genome
    * @param minCov     Minimum alignment coverage
    * @param id2seqname Map of id -> genome/chrm/contig/scaffold name
    * @return Syntenic anchors as 2-tuple of gene IDs
    */
  def findSytenicAnchors(fingertrees: Map[String, FingerTree],
                         id2seqname: Map[Int, (String, Int)],
                         minCov: Double,
                         splitOri: Boolean,
                         alignments: List[(Int, List[List[Alignment]])]
                        ): List[(Gene, Gene)] = {
    //set method to filter alignments by given minimum alignment coverage
    val alignment_filter = filterByCoverage(minCov)

    @tailrec def _findSyntenicAnchors(remaining_alignments: List[(Int, List[List[Alignment]])],
                                      syntenic_anchors: List[(Gene, Gene)]
                                     ): List[(Gene, Gene)] = {

      //check whether there are still alignments to process
      remaining_alignments match {
        //no alignments to process
        case Nil => syntenic_anchors
        //alignments to process
        case current_alignment :: tail => {
          //get sequence ID of the current query
          val query_seqname = id2seqname(current_alignment._1)._1
          //println("FINDING SA FOR: " + query_seqname)
          //get query's finger tree
          val query_finger_tree = fingertrees(query_seqname)
          //iterate through each alignment and construct ptolemy's gene graph
          val sa = current_alignment._2.foldLeft(syntenic_anchors)((acc_sa, alignment_interval) => {
            //get overlapping genes in current query's aligned sequence; using head since query coords are all the
            // same for each alignment interval
            val query_local_genes = alignment_filter(alignment_interval.head.qcoords,
              query_finger_tree.filterOverlaps(alignment_interval.head.qcoords).toList)
            //println("OVERLAPPING GENES: " + query_local_genes)
            //determine syntenic anchors by identifying overlapping genes between query and ref(s)
            query_local_genes.foldLeft(acc_sa)((interval_acc_sa, query_gene) => {
              //println("QUERYING GENE" + query_gene)
              //iterate through each alignment in the alignment interval
              alignment_interval.foldLeft(interval_acc_sa)((acc, alignment) => {
                //println("PROCESSING TARGET alignment: " + alignment)
                //project the coordinates of the current query gene onto the coordinates of the aligned ref sequence
                val projected_coords = (
                  normalizePosition(alignment.rcoords._1, alignment.rcoords._2, alignment.qcoords._1,
                    alignment.qcoords._2, query_gene._1._1),
                  normalizePosition(alignment.rcoords._1, alignment.rcoords._2, alignment.qcoords._1,
                    alignment.qcoords._2, query_gene._1._2)
                )
                //println("PROJECTED COORDS IN TARGET: " + projected_coords)
                //get overlapping genes with ref
                val overlapping_genes =
                alignment_filter(projected_coords, fingertrees(alignment.ref).filterOverlaps(projected_coords).toList)
                //println("TARGET OVERLAPPING GENES: " + overlapping_genes)
                //println(overlapping_genes)
                //no overlapping genes, move on
                if (overlapping_genes.isEmpty) acc
                //one or more overlapping gene
                else {
                  if (overlapping_genes.size > 1) {
                    println("Multiple genes forming syntenic anchor: projected coords for " +
                      "query gene " + query_gene + " from " + query_seqname + " is " + projected_coords + " and " +
                      "contains the following overlapping genes in " + alignment.ref + " " + overlapping_genes)
                  }
                  //user specified to collapse genes only when they have the same orientation, here they don't
                  if (splitOri && !overlapping_genes.forall(_._2.ori == query_gene._2.ori)) acc
                  //add syntenic anchors
                  else overlapping_genes.foldLeft(acc)((a, overlap) => (query_gene._2, overlap._2) :: a)
                }
              })
            })
          })
          //println(sa)
          //continue identifying syntenic anchors
          _findSyntenicAnchors(tail, sa)
        }
      }
    }

    _findSyntenicAnchors(alignments, List())
  }

  /**
    * Implementation of ptolemy.
    * Constructs gene graph based on whole genome alignment. There are 3 main steps:
    * 1). Finding syntenic anchors based gene projections from whole genome alignments
    * 2). Identifying ortholog clusters based on connected components from syntenic anchor graph
    * 3). Assignment of new canonical gene ID for each ortholog cluster
    *
    * @param fingertrees Map of genome -> fingertree data structure
    * @return 5-tuple as (GeneGraph, Paths, SA-mappings, Node coverage, Edge coverage). SA-mappings is map(old id ->
    *         new id) for the new assigned id for each ortholog cluster
    */
  def ptolemyGeneGraph(fingertrees: Map[String, FingerTree],
                       sa_graph: Map[Int, List[Int]],
                      ): (GeneGraph, Paths, Map[Int, Int], Map[Int, Int], Map[(Int, Int), Int]) = {
    println(timeStamp + "Identifying syntenic anchors")
    println(timeStamp + "--Identifying connected components")
    //fetch ortholog clusters, which are connected components
    val sa_ccs = getConnectedComponents(sa_graph)
    //get total number of ccs
    val sa_ccs_size = sa_ccs.size
    println(timeStamp + "--Found " + sa_ccs_size + " connected components")
    //define new canonical gene ID for each ortholog cluster/ccs
    val new_ids = defineCanonicalIDs(sa_ccs)
    //iterate through each fingertree and create gene graph
    val global_paths = {
      fingertrees.toList.foldLeft(empty_paths) { case (genome_paths, (genome, genes)) => {
        //update genome paths
        genome_paths + (
          //sort by gene id, get new gene ids, DO NOT collapse overlapping cases (for now)
          genome -> genes.toList.sortBy(_._2.id).map(x => new Gene(new_ids.getOrElse(x._2.id, x._2.id), '+')))
      }
      }.mapValues(_.toSet.toList)
    }
    //construct gene graph along with node and edge coverage
    val (global_gene_graph, node_coverage, edge_coverage) = {
      //iterate through each path and update accordingly
      val tmp = global_paths.foldLeft((empty_gene_graph, empty_node_coverage, empty_edge_coverage)){
        case ((graph, ncov, ecov), (genome, path)) => updateGeneGraph(path, false, graph, ncov, ecov)
      }
      (tmp._1.mapValues(_.distinct), tmp._2, tmp._3)
    }

    //total gene graph nodes
    val total_nodes = computeTotalNodes(global_gene_graph)

    println(timeStamp + "Constructed gene graph from whole-genome alignments of " + total_nodes + " nodes and " +
      global_gene_graph.toList.map(_._2.size).sum + " edges")
    //return global gene graph
    (global_gene_graph, global_paths, new_ids, node_coverage, edge_coverage)
  }
}
