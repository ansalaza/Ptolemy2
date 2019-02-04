package misc

import java.io.{File, PrintWriter}

import de.sciss.fingertree.RangedSeq
import utilities.SequenceFormatUtils.fetchSubSeqs
import utilities.FileHandling._
import utilities.IntervalUtils.isContained

/**
  * Author: Alex N. Salazar
  * Created on 23-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object FetchSubSequences {

  case class Config(
                     sequencesFile: File = null,
                     configFile: File = null,
                     prefix: String = null,
                     outputDir: File = null,
                     normalize: Boolean = false,
                     gffFile: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("fetch-sequences") {
      opt[File]('s', "sequences") required() action { (x, c) =>
        c.copy(sequencesFile = x)
      } text ("Reads/Assembly file in FASTQ/FASTA-formatted file.")
      opt[File]('c', "config-file") required() action { (x, c) =>
        c.copy(configFile = x)
      } text ("Tab-delimited file where each line contains: sequence ID and (optionally) start end coordinates. " +
        "Will fetch entire sequence if only the sequence ID is provided.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      note("\nOPTIONAL")
      opt[File]("gff") action { (x, c) =>
        c.copy(gffFile = x)
      } text ("Extract genes from given GFF file contained within sub-sequences.")
      opt[Unit]("normalize") action { (x, c) =>
        c.copy(normalize = true)
      } text ("Normalize coordinates of genes to reflect extracted sub-sequences. Assumes '--gff' parameter.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.sequencesFile)
      getSequences(config)
    }
  }

  def getSequences(config: Config): Unit = {
    println(timeStamp + "Processing configuration file")
    //load sequences to process
    val config_map = openFileWithIterator(config.configFile).foldLeft(List[(String, Option[(Int, Int)])]()) {
      (map, line) => {
        //split line
        val split = line.split("\t")
        //only sequence id was provided, add to list
        if (split.size == 1) (split.head, None) :: map
        //sequence id and coordinates were provided
        else {
          //sanity check
          assert(split.size >= 3, "Unexpected number of columns in the following line: " + line)
          //ad to list
          (split.head, Option((split(1).toInt, split(2).toInt))) :: map
        }
      }
    } //group by sequence identifier
      .groupBy(_._1)
      //create a list of all coordinates, if any, and sort by starting coordinate
      .mapValues(_.foldLeft(List[(Int, Int)]())((b, a) => if (a._2.isEmpty) b else a._2.get :: b).sortBy(_._1))
    println(timeStamp + "Found " + config_map.size + " (sub)sequences to extract")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".fa")
    //output (sub)sequences to disk
    fetchSubSeqs(config_map, pw, config.sequencesFile)
    pw.close

    if (config.gffFile != null) {
      //create a map of seq ID -> fingertree for all the extracted subseq coordinates
      val ft_map = {
        //set empty fingertree
        val empty = RangedSeq.empty[((Int, Int), Int), Int](_._1, Ordering.Int)
        //add subseq coords to fingertree
        config_map.foldLeft((Map[String, RangedSeq[((Int, Int), Int), Int]](), 0)) {
          //iterate through each coord in seq
          case ((map, id), (seq, coords)) => coords.foldLeft((map, id)) { case ((acc_map, acc_id), coord) => {
            //get current finger tree
            val current = acc_map.getOrElse(seq, empty)
            //add coord, exclusive end coord
            (acc_map + (seq -> (current + ((coord._1, coord._2 + 1), id))), acc_id + 1)
          }
          }
        }._1
      }
      //create output file
      val pw2 = new PrintWriter(config.outputDir + "/" + config.prefix + ".gff")
      println(timeStamp + "Parsing GFF file")
      openFileWithIterator(config.gffFile).foreach(line => {
        val columns = line.split("\t")
        val coords = ft_map.get(columns.head)
        // sequence exists in config map
        if (coords.nonEmpty) {
          //take all annotations in this sequence
          if (coords.get.isEmpty) pw2.println(line)
          //take only those within coords
          else {
            //get coords of current annotation
            val currrent_coords = (columns(3).toInt, columns(4).toInt)
            //get coords from config map that are contained
            val contained = coords.get.filterIncludes((currrent_coords)).toList
            //annnotation is contained
            if (contained.nonEmpty) {
              assert(contained.size == 1)
              if (!config.normalize) pw2.println(line)
              else pw2.println(columns.zipWithIndex.map { case (column, index) => {
                if(index == 0) column + "_" + contained.head._1._1 + "_" + (contained.head._1._2 - 1)
                else if (!(index == 3 || index == 4)) column
                else (column.toInt - contained.head._1._1).toString
              }
              }.mkString("\t"))
            }
          }
        }
      }
      )
      pw2.close
      println(timeStamp + "Successfully completed!")
    }
  }
}
