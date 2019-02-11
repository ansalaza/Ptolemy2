package utilities

import java.io.{File, PrintWriter}

import utilities.FileHandling.openFileWithIterator
import utilities.GFFutils.Gene
import utilities.GeneGraphUtils.{GeneGraph}

/**
  * Author: Alex N. Salazar
  * Created on 14-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GFAutils {

  /**
    * Set of methods for writing a file to GFA-format
    */
  trait GFAwriter {

    /**
      * Method to create a path line in GFA format
      *
      * @param path_name
      * @param path
      * @return
      */
    def makePathLine(path_name: String, path: List[Gene], length: Int): String = {
      val path_string = path.map(x => x.id.toString + x.ori).mkString(",")
      "P\t" + path_name + "\t" + path_string + "\t" + length
    }

    /**
      * Method to create segment line for ptolemy's gfa format
      * @param node
      * @param gcov
      * @param acov
      * @return
      */
    def makeSegmentLine(node: Int, gcov: Int, acov: Int): String = "S\t" + node + "\t\tGC:i:" + gcov + "\tRC:i:" + acov

    /**
      * Method to create link line for ptolemy's gfa format
      * @param edge
      * @param gcov
      * @param acov
      * @return
      */
    def makeLinkLine(edge: (Int, Int), gcov: Int, acov: Int): String = {
      "L\t" + edge._1 + "\t+\t" + edge._2 + "\t+\t\t" + "GC:i:" + gcov + "\tRC:i:" + acov
    }

    /**
      * Method to output ptolemy gene graph to GFA format (ptolemy version)
      * @param pw
      * @param graph
      * @param paths
      * @param genome_node_coverage
      * @param genome_edge_coverage
      * @param alignment_node_coverage
      * @param alignment_edge_coverage
      * @param tmp_paths
      */
    def ptolemyGraph2GFA(pw: PrintWriter,
                         graph: GeneGraph,
                         paths: Map[String, List[Gene]],
                         genome_node_coverage: Map[Int, Int],
                         genome_edge_coverage: Map[(Int,Int), Int],
                         alignment_node_coverage: Map[Int,Int],
                         alignment_edge_coverage: Map[(Int,Int), Int],
                         tmp_paths: File = null): Unit = {

      //output header
      pw.println("H\tGene-graph")
      //get all nodes and report cov
      graph.keys.toList.sorted.foreach(node => {
        //output line
        pw.println(
          makeSegmentLine(node,
            genome_node_coverage.getOrElse(node, 0),
            alignment_node_coverage.getOrElse(node, 0)))
      })
      //iterate through each node and report edge coverage
      graph.toList.foreach { case (node, edges) => {
        //iterate through each edge
        edges.foreach(edge => {
          //output line
          pw.println(
            makeLinkLine((node, edge),
              genome_edge_coverage.getOrElse((node, edge), 0),
              alignment_edge_coverage.getOrElse((node, edge), 0)
            )
          )
        })
      }}
      //iterate through each path and report
      paths.toList.foreach { case (genome, path) => pw.println(makePathLine(genome, path, -1)) }
      //add read paths, if they exist
      if(tmp_paths != null) openFileWithIterator(tmp_paths).foreach(line => pw.println(line))
    }

    /**
      * Output a given gene-graph to disk in GFA format
      *
      *
    def graph2GFA(pw: PrintWriter,
                  graph: GeneGraph,
                  paths: Map[String, List[Int]): Unit = {
      //get coverage for every node
      val node_coverage = paths.toList.foldLeft(Map[Int, Int]()) { case (cov_map, (genome, path)) => {
        //iterate through each path and update coverage
        path.foldLeft(cov_map)((acc_cov_map, node) => {
          //get coverage of current node
          val cov = acc_cov_map.getOrElse(node, 0) + 1
          //update coverage
          acc_cov_map + (node -> cov)
        })
      }
      }
      //get coverage for every edge
      val edge_coverage = paths.toList.foldLeft(Map[(Int, Int), Int]()) { case (cov_map, (genome, path)) => {
        path.sliding(2).foldLeft(cov_map)((acc_cov_map, _edge) => {
          //get forward edge
          val forward_edge = (_edge.head, _edge(1))
          //get reverse edge
          val reverse_edge = forward_edge.swap
          //get forward edge coverage
          val forward_edge_cov = acc_cov_map.getOrElse(forward_edge, 0) + 1
          //get reverse edge
          val reverse_edge_cov = acc_cov_map.getOrElse(reverse_edge, 0) + 1
          //update edge
          (acc_cov_map + (forward_edge -> forward_edge_cov)) + (reverse_edge -> reverse_edge_cov)
        })
      }
      }

      //output header
      pw.println("H\tGene-graph")
      //get all nodes and report cov
      paths.values.toList.flatten.toSet.toList.sorted
        .foreach(node => pw.println(makeSegmentLine(node, node_coverage(node))))
      //iterate through each node and report edge coverage
      graph.toList.foreach { case (node, edges) => {
        edges.foreach(edge => pw.println(makeLinkLine((node, edge), edge_coverage.getOrElse((node, edge),0))))
      }
      }
      paths.toList.foreach { case (genome, path) => pw.println(makePathLine(genome, path)) }
    }
      */
  }

  trait GFAreader {

    /**
      * Function to parse a path line from a GFA file.
      * @return 2-tuple as (Paths, optional size)
      */
    def parsePathLine: String => (String, List[Gene], Int) = line => {
      //split line
      val split = line.split("\t")
      assert(split.size >= 2, "Unexpected number of columns in line: " + line)
      //get path
      val path = {
        if(split(2).isEmpty) List[Gene]()
        else {
          split(2).split(",").toList.map(x => {
            val (node, ori) = x.partition(_.isDigit)
            assert(ori.size == 1 && (ori.head == '-' || ori.head == '+'), "Unexpected orientation for path line: " + line)
            new Gene(node.toInt,ori.head)
          })
        }
      }
      //get size, if available
      //val size = if(split.size >= 4) None else Option(split(3).toInt)
      //return path and optional size
      (split(1), path, split(3).toInt)
    }

    /**
      * Function to load GFA file given a file object
      * @return 2-tuple (GeneGraph and List of paths)
      */
    def loadGFA: File => (GeneGraph, Map[String, List[Gene]]) = file => {
      val (all_nodes, all_edges, all_paths) = {
        openFileWithIterator(file).foldLeft((List[Int](), List[(Int, Int)](), List[(String, List[Gene])]())) {
          case ((nodes, edges, paths), line) => {
            val columns = line.split("\t")
            columns.head match {
              case "S" => (columns(1).toInt :: nodes, edges, paths)
              case "L" => {
                val (node1, node2) = (columns(1).toInt, columns(3).toInt)
                (nodes, (node1, node2) :: edges, paths)
              }
              case "P" => {
                val (id, path, size) = parsePathLine(columns(2))
                (nodes, edges, (id, path) :: paths)
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
      //sanity check that there are bi-direction edges
      graph.toList.foreach{case (node, edges) => {
        edges.foreach(edge => {
          assert(graph(edge).toSet.contains(node), "Expected bi-directional edge for " + (node,edge))
        })
      }}
      (graph, all_paths.toMap)
    }


  }

  trait GFAconverter {
    /**
      * Method to convert GFA-formatted file to CSV file
      * @param gfa
      * @param pw_node
      * @param pw_edge
      */
    def gfa2CSV(gfa: File, pw_node: PrintWriter, pw_edge: PrintWriter): Unit = {
      //output headers
      pw_node.println("Source;Weight")
      pw_edge.println("Source;Target;Weight")
      //iterate through each line in gfa file
      openFileWithIterator(gfa).toList.foreach(line => {
        //split
        val split = line.split("\t")
        //attempt to get coverage
        val cov = split.find(_.startsWith("FC:i"))
        //only process link lines
        split.head match {
          //output segment line
          case "S" => pw_node.println(split(1) + (if (cov.isEmpty) "" else ";" + cov.get.filter(_.isDigit)))
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
