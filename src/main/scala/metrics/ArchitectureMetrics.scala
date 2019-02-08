package metrics

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.GFAutils.GFAreader
import utilities.GFFutils.Gene
import utilities.NumericalUtils.max
import atk.ProgressBar.progress
import scala.annotation.tailrec

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
                     minLength: Int = -1,
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
      note("\nOPTIONAL\n")
      opt[Int]("min-length") action { (x, c) =>
        c.copy(minLength = x)
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
    println(timeStamp + "Processing paths from GFA file")
    //fetch architectures from GFA file as 2-tuple: (architecture/path, (total instances, list of all names))
    val architectures = {
      //iterate through each line and process only path lines
      openFileWithIterator(config.gfaFile).foldLeft(List[(String, List[Gene])]())((paths, line) => {
        //split line
        val columns = line.split("\t")
        //check line type
        columns.head match {
          //path line
          case "P" => {
            //parse line and get path name and path
            val (name, path, size) = parsePathLine(line)
            //only add if the path is not empty or is large enough
            if (path.isEmpty || (config.minLength != -1 && size < config.minLength)) paths
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
        //group by architecture
        .groupBy(_._2)
        //sort by most frequent
        .toList.sortBy(- _._2.size)
        //assign unique IDs as keys and values as 4-tuples: (Path, Ori value, frequency, Seq IDs)
        .foldLeft((Map[Int, (List[Gene], Int, Int, List[String])](), 0)) { case ((map, id), path) => {
        (map + ((id) -> (path._1, path._1.map(_.ori).toSet.size - 1, path._2.size, path._2.map(_._1))), id + 1)
      }
      }._1
    }
    println(timeStamp + "Found a total of " + architectures.size + " unique paths")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    //output header
    pw.println("ID\tPath\tSingleOri\tCount\tNames")
    //set order for columns for metrics and distance matrix
    val rows = architectures.keys.toList.sorted
    //iterate through each architecture and output to disk
    rows.foreach(id => {
      val (path, ori, count, names) = architectures(id)
      pw.println(
        id + "\t" + path.map(x => x.id.toString + x.ori).mkString(",") + "\t" + ori + "\t" + count +
          "\t" + names.mkString(","))
    })
    pw.close()

    /**
      * Function to compute the normalized edit distance of two given architecture IDs. This is defined as:
      * dist = min{ editDist(x,y), editDist(x, y.reverse), editDist(x.reverse, y), editDist(x.reverse, y.reverse }
      * size = max{ x.size, y.size }
      * 1 - ( dist / size)
      *
      * @return Double
      */
    def computeNormEditDist: (Int, Int) => Double = (_x, _y) => {
      //get architectures
      val (x, y) = (architectures(_x)._1.map(_.id), architectures(_y)._1.map(_.id))
      //get max size
      val max_size = max(x.size, y.size).toDouble
      //no overlap, return max normalized distance
      if(x.toSet.intersect(y.toSet).size == 0) 1.0
      //compute distance for all possible forms
      else {
        //compute edit distance as:
        val edit_distance = List(
          //normal orientations
          editDist(x, y),
          //y-reverse
          editDist(x, y.reverse),
          //x-reverse
          editDist(x.reverse, y),
          //both-reverse
          editDist(x.reverse, y.reverse)
        ).min
        //return normalized edit distance
        (edit_distance / max_size)
      }
    }

    /**
      * Create a map of all pairwise distances (upper-diaganol only)
      */
    val dist_map = {
      val all_pairwise = computeAllPairwise(architectures.keys.toList, List())
      println(timeStamp + "Computing " + all_pairwise.size + " pairwise path-distances")
      all_pairwise.map { case (subj, target) => {
        progress(100000)
        ((subj, target), computeNormEditDist(subj, target))
      }}.toMap
    }
    //set output file
    val pw2 = new PrintWriter(config.outputDir + "/" + config.prefix + ".distances.matrix")
    pw2.println("$" + "\t" + rows.mkString("\t"))
    //create distance matrix
    rows.foreach(row => {
      //output row name
      pw2.print(row)
      rows.foreach(column => {
        if(row == column) pw2.print("\t" + 0.0)
        else pw2.print("\t" + dist_map.getOrElse((row, column), dist_map((column, row))))
      })
      pw2.println
    })
    pw2.close
    println(timeStamp + "Successfully completed!")
  }

  /**
    * Function to compare a given path and it's reverse orientation and return the smallest of the to
    *
    * @return List[Gene]
    */
  def getSmallestPath: List[Gene] => List[Gene] = forward => {
    //get forward and reverse orientation
    val reverse = forward.reverse.map(_.reverse())
    //set total path size
    val total_size = forward.size

    /**
      * Get smallest path of forward and reverse orientation
      *
      * @param i Index
      * @return List[Gene]
      */
    def _getSmallestPath(i: Int): List[Gene] = {
      if (i == total_size) forward
      else if (forward(i).id < reverse(i).id) forward
      else if (forward(i).id > reverse(i).id) reverse
      else _getSmallestPath(i + 1)
    }

    //get smallest path
    _getSmallestPath(0)
  }

  /**
    * Compute all possible pairwise interactions (upper-diagonal only) given a list of IDs
    *
    * @param ids
    * @param interactions
    * @return
    */
  @tailrec def computeAllPairwise[A](ids: List[A], interactions: List[(A, A)]): List[(A, A)] = {
    ids match {
      //no more pairwise interactions to get
      case Nil => interactions
      //remainin interactions
      case (subj :: tail) => {
        //iterate and obtain upper diagonal for current id
        computeAllPairwise(tail,
          tail.foldLeft(interactions)((acc, target) => if (subj == target) acc else (subj, target) :: acc),
        )
      }
    }
  }

  /**
    * Edit distance implementation. Thanks! https://www.reddit.com/r/scala/comments/7sqtyf/scala_edit_distance_implementation/
    *
    * @param a
    * @param b
    * @tparam A
    * @return
    */
  def editDist[A](a: Iterable[A], b: Iterable[A]): Int = {
    val startRow = (0 to b.size).toList
    a.foldLeft(startRow) { (prevRow, aElem) =>
      (prevRow.zip(prevRow.tail).zip(b)).scanLeft(prevRow.head + 1) {
        case (left, ((diag, up), bElem)) => {
          val aGapScore = up + 1
          val bGapScore = left + 1
          val matchScore = diag + (if (aElem == bElem) 0 else 1)
          List(aGapScore, bGapScore, matchScore).min
        }
      }
    }.last
  }
}
