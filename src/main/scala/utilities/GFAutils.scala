package utilities

import java.io.{File, PrintWriter}

import utilities.FileHandling.openFileWithIterator
import utilities.GFFutils.Gene
import utilities.GeneGraphUtils.{GeneGraph, GenePaths}
import utilities.SequenceGraphUtils.{SequenceGraph, SequencePaths}

/**
  * Author: Alex N. Salazar
  * Created on 14-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GFAutils {

  /**
    * Method to create a path line in GFA format
    *
    * @param path_name
    * @param path
    * @return
    */
  def makePathLine(path_name: String, path: List[String], length: Int): String = {
    "P\t" + path_name + "\t" + path.mkString(",") + "\t" + length
  }


  /**
    * Method to create segment line for ptolemy's gfa format
    *
    * @param node
    * @param cov
    * @return
    */
  def makeSegmentLine(node: Any, seq: Option[String], cov: Int): String =
    "S\t" + node + "\t" + (if (seq.isEmpty) "" else seq.get) + "\tFC:i:" + +cov

  /**
    * Method to create link line for ptolemy's gfa format
    *
    * @param edge
    * @param cov
    * @return
    */
  def makeLinkLine(edge: (Any, Any), cov: Int): String = "L\t" + edge._1 + "\t+\t" + edge._2 + "\t+\t\tFC:i:" + cov

  trait SequenceGraphWriter {
    def seqGraph2GFA(pw: PrintWriter,
                     graph: SequenceGraph,
                     paths: SequencePaths,
                     node_coverage: Map[String, Int],
                     edge_coverage: Map[(String, String), Int]): Unit = {

      //output header
      pw.println("H\tSequence-graph")
      //get all nodes and report cov
      graph.keys.toList.sortBy(x => (x.seq_id, x.gene_id)).foreach(node => {
        pw.println(makeSegmentLine(node.gfaID(), Option(node.aa.toString), node_coverage.getOrElse(node.gfaID(), 0)))
      })

      //iterate through each node and report edge coverage
      graph.toList.foreach { case (node, edges) => {
        //iterate through each edge
        edges.foreach(edge => {
          //output line
          pw.println(makeLinkLine((node.gfaID(), edge.gfaID()),
            edge_coverage.getOrElse((node.gfaID(), edge.gfaID()), 0)))
        })
      }
      }
      //iterate through each path and report
      paths.toList.foreach { case (genome, path) =>
        pw.println(makePathLine(genome, path.map(x => x.gfaID() + "+"), path.size))
      }
    }

  }


  /**
    * Set of methods for writing a file to GFA-format
    */
  trait GeneGraphWriter {


    /**
      * Method to output ptolemy gene graph to GFA format (ptolemy version)
      *
      * @param pw
      * @param graph
      * @param alignment_node_coverage
      * @param alignment_edge_coverage
      * @param tmp_paths
      */
    def geneGraph2GFA(pw: PrintWriter,
                      graph: GeneGraph,
                      paths: Option[GenePaths],
                      alignment_node_coverage: Map[Int, Int],
                      alignment_edge_coverage: Map[(Int, Int), Int],
                      tmp_paths: File = null): Unit = {

      //output header
      pw.println("H\tGene-graph")
      //get all nodes and report cov
      graph.keys.toList.sorted.foreach(node => {
        //output line
        pw.println(makeSegmentLine(node, None, alignment_node_coverage.getOrElse(node, 0)))
      })
      //iterate through each node and report edge coverage
      graph.toList.foreach { case (node, edges) => {
        //iterate through each edge
        edges.foreach(edge => {
          //output line
          pw.println(makeLinkLine((node, edge), alignment_edge_coverage.getOrElse((node, edge), 0))
          )
        })
      }
      }
      if(paths.nonEmpty) {
        paths.get.foreach(path => pw.println(makePathLine(path._1, path._2._1.map(_.gfaID()), path._2._2)))
      }
      //add read paths, if they exist
      if (tmp_paths != null) openFileWithIterator(tmp_paths).foreach(line => pw.println(line))
    }

  }

  trait GeneGraphReader {

    /**
      * Funtion to parse path column into a list of genes
      *
      * @return
      */
    def parsePath: String => List[Gene] = line => {
      line.split(",").toList.map(x => {
        val (node, ori) = x.partition(_.isDigit)
        assert(ori.size == 1 && (ori.head == '-' || ori.head == '+'), "Unexpected orientation for path line: " + line)
        new Gene(node.toInt, ori.head)
      })
    }

    /**
      * Function extract edge labels from a given Paths object
      *
      * @return Map[(Gene,Gene), List[String]
      */
    def paths2EdgeLabels: GenePaths => Map[(Gene, Gene), Set[String]] = paths => {
      paths.foldLeft(Map[(Gene, Gene), List[String]]()) { case (quiver, (name, (path, size))) => {
        path.size match {
          case 1 => quiver
          case _ => path.sliding(2).map(x => (x(0), x(1))).foldLeft(quiver)((acc, edge) => {
            //get forward labels
            val forward = acc.getOrElse(edge, List[String]())
            //get reverse labels
            val reverse = acc.getOrElse(edge.swap, List[String]())
            //update quiver
            (acc + (edge -> (name :: forward))) + (edge.swap -> (name :: reverse))
          })
        }
      }
      }.mapValues(_.toSet)
    }

    /**
      * Function to load GFA file given a file object
      *
      * @return 2-tuple (GeneGraph and List of paths)
      */
    def loadGFA: File => (GeneGraph, GenePaths) = file => {
      val (all_nodes, all_edges, all_paths) = {
        openFileWithIterator(file).foldLeft((List[Int](), List[(Int, Int)](), List[(String, (List[Gene], Int))]())) {
          case ((nodes, edges, paths), line) => {
            val columns = line.split("\t")
            columns.head match {
              case "S" => (columns(1).toInt :: nodes, edges, paths)
              case "L" => {
                val (node1, node2) = (columns(1).toInt, columns(3).toInt)
                (nodes, (node1, node2) :: edges, paths)
              }
              case "P" => {
                assert(columns.size == 4, "Unexpected number of columns in line: " + line)
                val (id, size) = (columns(1), columns.last.toInt)
                val path = parsePath(columns(2))
                (nodes, edges, (id, (path, size)) :: paths)
              }
              case _ => (nodes, edges, paths)
            }
          }
        }
      }
      //iterate through each edge and add to graph initially created with all nodes and empty edges
      val graph = all_edges.foldLeft(all_nodes.map(node => (node, List[Int]())).toMap)((graph, edge) => {
        //sanity check for source node
        assert(graph.get(edge._1).nonEmpty, "Unexpected node when appending edge for node: " + edge._1)
        //update graph
        (graph + (edge._1 -> (edge._2 :: graph(edge._1))))
      })
      //sanity check that all nodes are accounted for
      all_nodes.foreach(node => assert(graph.contains(node), "Expected node " + node + " in graph"))
      //sanity check that there are bi-direction edges
      graph.toList.foreach { case (node, edges) => {
        edges.foreach(edge => {
          assert(graph(edge).toSet.contains(node), "Expected bi-directional edge for " + (node, edge))
        })
      }
      }
      (graph, all_paths.toMap)
    }


  }

  trait GeneGraphConverter {
    /**
      * Method to convert GFA-formatted file to CSV file
      *
      * @param gfa
      * @param pw_node
      * @param pw_edge
      */
    def gfa2CSV(gfa: File,
                pw_node: PrintWriter,
                pw_edge: PrintWriter,
                labels: Map[Int, String] = Map()): Unit = {
      //output headers
      pw_node.println("Id;Weight" + (if (labels.isEmpty) "" else ";Label"))
      pw_edge.println("Source;Target;Weight")
      //iterate through each line in gfa file
      openFileWithIterator(gfa).toList.foreach(line => {
        //split
        val split = line.split("\t")
        //attempt to get coverage
        val cov = split.find(_.startsWith("FC:i:"))
        //only process link lines
        split.head match {
          //output segment line
          case "S" => pw_node.println(split(1) + (if (cov.isEmpty) "" else ";" + cov.get.filter(_.isDigit)) +
            (if (labels.isEmpty) "" else (";\"" + labels.getOrElse(split(1).toInt, "None")) + "\""))
          //output link line
          case "L" => {
            //output information in CSV-format
            pw_edge.println(split(1) + ";" + split(3) + (if (cov.isEmpty) "" else ";" + cov.get.filter(_.isDigit)))
          }
          case _ => "nothing"
        }
      })
    }
  }

}
