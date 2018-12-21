package utilities

import java.io.{File, PrintWriter}

import utilities.FileHandling.openFileWithIterator
import utilities.GeneGraphUtils.GeneGraph

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
      * @param _length
      * @return
      */
    def makePathLine(path_name: String, path: List[Int], _length: Option[Int] = None): String = {
      val length = if (_length.isEmpty) "" else _length.get
      "P\t" + path_name + "\t" + path.mkString(",") + "\t" + length
    }

    /**
      * Method to make segment line in GFA format
      *
      * @param node
      * @param cov
      * @return
      */
    def makeSegmentLine(node: Int, cov: Int): String = "S\t" + node + "\t\tFC:i:" + cov

    /**
      * Method to make link line in GFA format
      *
      * @param edge
      * @param cov
      * @return
      */
    def makeLinkLine(edge: (Int, Int), cov: Int): String = {
      "L\t" + edge._1 + "\t+\t" + edge._2 + "\t+\t\t" + "FC:i:" + cov
    }

    /**
      * Output a given gene-graph to disk in GFA format
      *
      * @param pw
      * @param graph
      * @param paths
      */
    def graph2GFA(pw: PrintWriter, graph: GeneGraph, paths: Map[String, List[Int]]): Unit = {
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
  }

  trait GFAreader {

    /**
      * Function to parse a path line from a GFA file.
      * @return 2-tuple as (Paths, optional size)
      */
    def parsePathLine: String => (List[Int], Option[Int]) = line => {
      //split line
      val split = line.split("\t")
      assert(split.size >= 2, "Unexpected number of columns in line: " + line)
      //get path
      val path = split(2).split(",").toList.map(_.toInt)
      //get size, if available
      val size = if(split.size >= 4) None else Option(split(3).toInt)
      //return path and optional size
      (path, size)
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
