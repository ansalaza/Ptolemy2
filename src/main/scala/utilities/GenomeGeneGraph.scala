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
    val (global_gene_graph, global_genome_paths) = {
      val tmp = fingertrees.toList.foldLeft((empty_gene_graph, empty_paths)) {
        case ((gene_graph, genome_paths), (genome, genes)) => {
          //get current genome path using new gene ids
          val local_genome_paths = {
            //sort by gene id, get new gene ids, DO NOT collapse overlapping cases (for now)
            genes.toList.sortBy(_._2.id).map(x => new Gene(new_ids.getOrElse(x._2.id, x._2.id), '+'))
          }
          //get all genes for current genome, get all edges and update gene graph
          val local_gene_graph = {
            if (local_genome_paths.size == 1) gene_graph + (local_genome_paths.head.id -> List())
            else local_genome_paths.sliding(2).foldLeft(gene_graph)((acc_gene_graph, _edge) => {
              //get new node IDs for the current two nodes
              val (node1, node2) = (_edge.head.id, _edge(1).id)
              //update node1's current edges
              val node1_current_edges = node2 :: acc_gene_graph.getOrElse(node1, List[Int]())
              //update node2's current edges
              val node2_current_edges = node1 :: acc_gene_graph.getOrElse(node2, List[Int]())
              //update map
              (acc_gene_graph + (node1 -> node1_current_edges)) + (node2 -> node2_current_edges)
            })
          }
          //return updated gene graph and update genome paths
          (local_gene_graph, genome_paths + (genome -> local_genome_paths))
        }
      }
      (tmp._1.mapValues(_.toSet.toList), tmp._2)
    }
    //total gene graph nodes
    val total_nodes = {
      //iterate through each node and it's corresponding edges
      global_gene_graph.foldLeft(List[Int]()) { case (all_nodes, (node, edges)) => {
        //add node to edges, iterate and add to all nodes
        (node :: edges).foldLeft(all_nodes)((b, a) => a :: b)
      }
      }.toSet.size
    }
    //compute node and edge coverage
    val (node_coverage, edge_coverage) = {
      //iterate through each genome
      global_genome_paths.toList.foldLeft((Map[Int, Int](), Map[(Int, Int), Int]())) {
        case ((node_cov, edge_cov), (genome, path)) => {
          //update node coverage
          val updated_node_coverage =
            path.foldLeft(node_cov)((cov, gene) => cov + (gene.id -> (cov.getOrElse(gene.id, 0) + 1)))
          //update edge coverage, both forward and reverse
          val updated_edge_coverage = {
            if (path.size == 1) edge_cov
            //iterate through each edge
            else path.sliding(2).foldLeft(edge_cov)((cov, _edge) => {
              //set forward edge
              val forward = (_edge.head.id, _edge(1).id)
              //self edge
              if (forward._1 == forward._2) {
                //updated forward coverage +2
                val updated_cov = cov.getOrElse(forward, 0) + 2
                cov + (forward -> updated_cov)
              } else {
                //updated forward coverage
                val updated_forward_cov = cov.getOrElse(forward, 0) + 1
                //set reverse edge
                val reverse = forward.swap
                //updated reverse edge
                val updated_reverse_cov = cov.getOrElse(reverse, 0) + 1
                //update coverage
                (cov + (forward -> updated_forward_cov)) + (reverse -> updated_reverse_cov)
              }
            })
          }
          (updated_node_coverage, updated_edge_coverage)
        }
      }
    }

    println(timeStamp + "Constructed gene graph from whole-genome alignments of " + total_nodes + " nodes and " +
      global_gene_graph.toList.map(_._2.size).sum + " edges")
    //return global gene graph
    (global_gene_graph, global_genome_paths, new_ids, node_coverage, edge_coverage)
  }
}
