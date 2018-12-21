package gene_graph

import java.io.{File, PrintWriter}
import utilities.FileHandling.{verifyFile, verifyDirectory, getFileName}
import utilities.GFAutils.GFAconverter

/**
  * Author: Alex N. Salazar
  * Created on 11-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GFAconverter extends GFAconverter {

  case class Config(
                     gfaFile: File = null,
                     outputDir: File = null,
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
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.gfaFile)
      gfa2Format(config)
    }
  }

  def gfa2Format(config: Config): Unit = {
    //parse formats
    config.formats.split("\\s+").foreach(format => {
      format.toLowerCase match {
        case "csv" => {
          //create output file
          val pw_node = new PrintWriter(config.outputDir + "/" + getFileName(config.gfaFile) + ".nodes.csv")
          val pw_edge = new PrintWriter(config.outputDir + "/" + getFileName(config.gfaFile) + ".edges.csv")
          //convert to CSV format
          gfa2CSV(config.gfaFile, pw_node, pw_edge)
          //close files
          pw_edge.close
          pw_node.close
        }
        case _ => println("WARNING: Unrecognize file format: " + format)
      }
    })
  }

}


