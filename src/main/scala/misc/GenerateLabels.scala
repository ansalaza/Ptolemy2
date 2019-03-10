package misc

import java.io.File

import utilities.FileHandling.{verifyDirectory, verifyFile, openFileWithIterator}

/**
  * Author: Alex N. Salazar
  * Created on 9-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GenerateLabels {

  case class Config(
                     graph: File = null,
                     outputDir: File = null,
                     path: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("generate-labels") {
      opt[File]('g', "graph") required() action { (x, c) =>
        c.copy(graph = x)
      } text ("Gene-graph in GFA-format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      opt[String]("path") required() action { (x, c) =>
        c.copy(path = x)
      } text ("Comma-separated string of all path IDs to generate labels for.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.graph)
      generateLabels(config)
    }
  }

  def generateLabels(config: Config): Unit = {
    //get path IDs
    val paths = config.path.replaceAll("\\s+","").split(",").toSet
    println("Source;Target;Type")
    //iterate through gfa and only process paths
    openFileWithIterator(config.graph).foreach(line => {
      val columns = line.split("\t")
      columns.head match {
        case "P" => {
          if(paths(columns(1))){
            val ids = columns(2).split("[+|-],")
            if(ids.size > 1) ids.sliding(2).foreach(x => println(List(x(0),x(1),columns(1)).mkString(";")))
          }
        }
        case _ => ()
      }
    })
  }

}
