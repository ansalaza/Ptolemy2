package misc

import java.io.{File, PrintWriter}
import utilities.SequenceFormatUtils.fetchSubSeqs
import utilities.FileHandling._

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
                     outputDir: File = null)

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
      note("OPTIONAL")
      opt[String]("prefix") action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
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
    }//group by sequence identifier
      .groupBy(_._1)
      //create a list of all coordinates, if any, and sort by starting coordinate
      .mapValues(_.foldLeft(List[(Int,Int)]())((b,a) => if(a._2.isEmpty) b else a._2.get :: b).sortBy(_._1))
    println(timeStamp + "Found " + config_map.size + " (sub)sequences to extract")
    //create output file
    val pw = {
      //get output name
      val name = if (config.prefix == null) getFileName(config.configFile) + ".subsequences" else config.prefix
      //sanity check for not overwriting an existing file
      assert(!new File(config.outputDir + "/" + name + ".fa").exists(), "Output file already exist: " + name)
      new PrintWriter(config.outputDir + "/" + name + ".fa")
    }
    //output (sub)sequences to disk
    fetchSubSeqs(config_map, pw, config.sequencesFile)
    pw.close
    println(timeStamp + "Successfully completed!")
  }
}
