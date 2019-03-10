package misc

import java.io.{File, PrintWriter}

import utilities.FileHandling.{getFileName, verifyDirectory, verifyFile, openFileWithIterator}
import utilities.GFAutils.GeneGraphConverter

/**
  * Author: Alex N. Salazar
  * Created on 11-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GFAconverter extends GeneGraphConverter {

  case class Config(
                     gfaFile: File = null,
                     outputDir: File = null,
                     coverage: Char = 'a',
                     labels: File = null,
                     formats: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("gfa-converter") {
      opt[File]('g', "gfa-file") required() action { (x, c) =>
        c.copy(gfaFile = x)
      } text ("Gene graph in GFA format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      opt[String]('f', "formats") required() action { (x, c) =>
        c.copy(formats = x)
      } text ("Comma-separated string containing format IDs to convert GFA to. Current supported options: 'csv'.")
      opt[Char]("coverage") required() action { (x, c) =>
        c.copy(coverage = x)
      } text ("Use [g]enome coverage or [r]ead coverage as node/edge weights.")
      note("\nOPTIONAL")
      opt[File]("labels") action { (x, c) =>
        c.copy(labels = x)
      } text ("Tab-delimited file containing: node ID, label.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.gfaFile)
      assert(config.coverage == 'g' || config.coverage == 'r', "Specify 'g' for genome coverage or 'r' for read " +
        "coverage.")
      gfa2Format(config)
    }
  }

  def gfa2Format(config: Config): Unit = {
    //load labels, if any
    val labels = {
      if(config.labels == null) Map.empty[Int,String]
      else openFileWithIterator(config.labels).toList.map(x =>{
        val c = x.split("\t")
        c(0).toInt -> c(1)
      }).toMap
    }
    //parse formats
    config.formats.split("\\s+").foreach(format => {
      format.toLowerCase match {
        case "csv" => {
          //create output file
          val pw_node = new PrintWriter(config.outputDir + "/" + getFileName(config.gfaFile) + ".nodes.csv")
          val pw_edge = new PrintWriter(config.outputDir + "/" + getFileName(config.gfaFile) + ".edges.csv")
          //convert to CSV format
          gfa2CSV(config.gfaFile, (if(config.coverage == 'g') "GC:i:" else "RC:i:"), pw_node, pw_edge, labels)
          //close files
          pw_edge.close
          pw_node.close
        }
        case _ => println("WARNING: Unrecognize file format: " + format)
      }
    })
  }

}


