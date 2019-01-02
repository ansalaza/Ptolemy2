package metrics

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.GFAutils.GFAreader
import utilities.GeneGraphUtils.PathEntry

/**
  * Author: Alex N. Salazar
  * Created on 27-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ArchitectureMetrics extends GFAreader {

  case class Config(
                     gfaFile: File = null,
                     outputDir: File = null,
                     prefix: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("path-metrics") {
      opt[File]('g', "gfa-file") required() action { (x, c) =>
        c.copy(gfaFile = x)
      } text ("Gene graph in GFA format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether input file/directories exists
      verifyDirectory(config.outputDir)
      verifyFile(config.gfaFile)
      geneGraphMetrics(config)
    }
  }

  def geneGraphMetrics(config: Config): Unit = {
    /**
      * Compare a given path and it's reverse form and return the smallest of the two
      * @return List[PathEntry]
      */
    def getSmallestPath: List[PathEntry] => List[PathEntry] = path => {
      //get hashcodes for forward and reverse
      val (forward,reverse) = (path, path.reverse.map(_.reverse()))
      //forward is smallest, returh path as is, else reverse it
      if(forward.hashCode() < reverse.hashCode()) forward else reverse
    }

    println(timeStamp + "Processing paths from GFA file")
    //fetch architectures from GFA file as 2-tuple: (architecture/path, (total instances, list of all names))
    val architectures = {
      //iterate through each line and process only path lines
      openFileWithIterator(config.gfaFile).foldLeft(List[(String, List[PathEntry])]())((paths, line) => {
        //split line
        val columns = line.split("\t")
        //check line type
        columns.head match {
          //path line
          case "P" => {
            //parse line and get path name and path
            val (name, path) = (columns(1), parsePathLine(line)._1)
            //only add if the path is not empty
            if (path.isEmpty) paths
            else {
              //get all different orientations of current path
              val path_orientations = path.map(_.ori).toSet.toList
              //check number of orientations
              val updated_path = path_orientations.size match {
                //only one orientation, reverse only if its in reverse orientation
                case 1 => if (path_orientations.head == '+') path else path.reverse.map(_.reverse())
                //dual orienation, retain smallest path
                case 2 => getSmallestPath(path)
              }
              (name, updated_path) :: paths
            }
          }
          //ignore all other lines
          case _ => paths
        }
      })
        //group by architecture, map values as 2-tuples (total paths, path names), sort by total paths
        .groupBy(_._2).mapValues(x => (x.size, x.map(_._1))).toList.sortBy(-_._2._1)
    }
    println(timeStamp + "Found a total of " + architectures.size + " unique paths")
    println(timeStamp + "Writing to disk")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    //output header
    pw.println("Path\tCount\tSingleOri\tNames")
    //iterate through each architecture and output to disk
    architectures.foreach { case (path, (count, names)) => {
      pw.println(
        path.map(x => x.nodeID.toString + x.ori).mkString(",") + "\t" +
        count + "\t" +
        (path.map(_.ori).toSet.size - 1) + "\t" +
        names.mkString(","))
    }
    }
    pw.close()
    println(timeStamp + "Successfully completed!")
  }


}
