package utilities

import java.io.{File, PrintWriter}

import utilities.AlignmentUtils.Alignment
import utilities.FileHandling.openFileWithIterator
import utilities.GFFutils.Gene
import utilities.GeneGraphUtils._
import utilities.GeneProjectionUtils.{FingerTree, filterByCoverage, projectGenes, adjustProjectionOri}

/**
  * Author: Alex N. Salazar
  * Created on 21-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object AlignmentGeneGraph {

  /**
    * Method to creage gene graph based on long-read alignments. Pseudocode:
    *
    * Initialize empty list of paths
    * Iterate for each read:
    * Initialize empty path for current read
    * For each alignment is sorted order:
    * Project genes onto alignment and add to path
    * Add path to list of paths
    * Initialize empty gene graph, node coverage, edge coverage
    * Iterate for each path:
    * update gene graph, node coverage, edge coverage
    *
    * @param fingertrees        Finger tree of reference genomes
    * @param alignments         List of all alignments per read
    * @param minCov             Minimum alignment coverage threshold
    * @param id2seqname         Map of ID -> (read name, size)
    * @param initial_gene_graph Initialized starting gene graph
    * @param sa_mapping         Syntenic anchor ID map containing old gene ID -> new gene ID
    * @param tmp_file           File to store temporary paths
    * @return 3-tuple (GeneGraph, Node coverage, Edge coverage, temp file of paths)
    */
  def alignmentGeneGraph(fingertrees: Map[String, FingerTree],
                         alignments: List[(Int, List[List[Alignment]])],
                         minCov: Double,
                         id2seqname: Map[Int, (String, Int)],
                         initial_gene_graph: GeneGraph,
                         sa_mapping: Map[Int, Int],
                         tmp_file: File
                        ): (GeneGraph, Map[Int, Int], Map[(Int, Int), Int]) = {
    //set method to filter alignments by given minimum alignment coverage
    val alignment_filter = filterByCoverage(minCov)
    //set curried method to filter out reads not covered by minimum alignment coverage
    val gene_projector = projectGenes(alignment_filter) _
    //set temporary projections file
    val pw = new PrintWriter(tmp_file)
    //iterate through each alignment interval and project genes on sequence
    alignments.foreach { case (read, alignment_interval) => {
      //obtain architecture on sequence
      val architecture = {
        //no alignments, return empty
        if (alignments.isEmpty) List[Gene]()
        //at least one alignment
        else {
          //iterate through each alignment interval and project genes
          val seq_projections = alignment_interval.map(alignments => {
            //project onto all alignments
            val projections = alignments.map(alignment =>
              (gene_projector(alignment, fingertrees(alignment.ref)), alignment.ori))
            //get largest projection
            //TODO: verify all projections with WGS alignments? what about partially-related overlapping seqs?
            //println(projections)
            //println(alignments)
            //println
            projections.sortBy(x => (x._1.size, if(x._1.isEmpty) Int.MaxValue else x._1.map(_.id).min)).head
          }).filter(_._1.nonEmpty)
          //no projections to make
          if(seq_projections.isEmpty) List[Gene]()
          //adjust orientation of final projection and replace gene id by new id assigned to respective ortholog cluster
          else adjustProjectionOri(seq_projections).map(x => new Gene(sa_mapping.getOrElse(x.id, x.id), x.ori))
        }
      }
      //get name and size
      val (name, size) = id2seqname(read)
      //create and output path line
      pw.println(makePathLine(name, architecture, size))
    }
    }
    //close tmp file
    pw.close()
    //iterate through each path and update gene graph as well as calculate node and edge coverage
    val (g, n, e) = {
      openFileWithIterator(tmp_file).foldLeft((empty_gene_graph, empty_node_coverage, empty_edge_coverage)) {
        case ((gene_graph, node_coverage, edge_coverage), line) => {
          //parse path line
          val (id, path, size) = parsePathLine(line)
          //update node coverage
          val updated_node_coverage =
            path.foldLeft(node_coverage)((cov, gene) => cov + (gene.id -> (cov.getOrElse(gene.id, 0) + 1)))
          //update edge coverage
          val (updated_gene_graph, updated_edge_coverage) = {
            //skip if only one gene
            if (path.size == 1) (gene_graph, edge_coverage)
            //iterate through path as edges, update coverage
            else path.sliding(2).foldLeft((gene_graph, edge_coverage)) { case ((graph, cov), _edge) => {
              //set nodes
              val (node1, node2) = (_edge.head.id, _edge(1).id)
              //set edge
              val edge = (node1, node2)
              //add node2 to current edges of node1
              (graph + (node1 -> (node2 :: graph.getOrElse(node1, List[Int]()))),
                //update coverage of edge
                cov + (edge -> (cov.getOrElse(edge, 0) + 1)))
            }
            }
          }
          //return updated graphs
          (updated_gene_graph, updated_node_coverage, updated_edge_coverage)
        }
      }
    }
    (g.mapValues(_.toSet.toList), n, e)
  }

}
