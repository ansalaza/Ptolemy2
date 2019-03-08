package utilities

import java.io.{File, PrintWriter}

import utilities.AlignmentUtils.Alignment
import utilities.FileHandling.openFileWithIterator
import utilities.GFFutils.Gene
import utilities.GeneGraphUtils._
import utilities.GeneProjectionUtils.{FingerTree, filterByCoverage, projectGenes}

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
            //println(alignments)
            //project onto all alignments
            val projections = alignments.map(alignment =>
              (gene_projector(alignment, fingertrees(alignment.ref)), alignment.ori))
            //println(projections.map(_._1.map(_.id).map(y => sa_mapping.getOrElse(y, y))))
            //get largest projection
            //TODO: verify all projections with WGS alignments? what about partially-related overlapping seqs?
            if (alignments.size > 1) {
              //println(projections)
              //println(alignments)
              //println
              (List[Gene](), '-')
            }
            else projections.head
          }).filter(_._1.nonEmpty)
          //no projections to make
          if (seq_projections.isEmpty) List[Gene]()
          //adjust orientation of final projection and replace gene id by new id assigned to respective ortholog cluster
          else seq_projections.flatMap(x => if (x._2 == '+') x._1 else x._1.reverse.map(_.reverse()))
            .map(x => new Gene(sa_mapping.getOrElse(x.id, x.id), x.ori))
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
      //load as gfa, but retain only paths
      loadGFA(tmp_file)._2.foldLeft((empty_gene_graph, empty_node_coverage, empty_edge_coverage)) {
        case ((gene_graph, node_coverage, edge_coverage), (name, (path,size))) => {
          //update with current path
          updateGeneGraph(path, true, gene_graph, node_coverage, edge_coverage)
        }
      }
    }
    (g.mapValues(_.distinct), n, e)
  }

}
