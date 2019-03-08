package misc

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.MetaDataUtils.loadDescriptions

/**
  * Author: Alex N. Salazar
  * Created on 20-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object FetchWithStrings {

  case class Config(
                     descriptionFile: File = null,
                     stringsFile: File = null,
                     outputDir: File = null,
                     geneMappingFile: File = null,
                     prefix: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("extract-genes") {
      opt[File]('d', "descriptions") required() action { (x, c) =>
        c.copy(descriptionFile = x)
      } text ("Descriptions file (see README).")
      opt[File]('s', "strings") required() action { (x, c) =>
        c.copy(stringsFile = x)
      } text ("File containing strings (case-insensitive) to identify in gene description, one per line. Comma " +
        "separated strings per line imply AND operator.")
      opt[String]('p', "prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output files.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      note("\nOPTIONAL")
      opt[File]("gene-families") action { (x, c) =>
        c.copy(geneMappingFile= x)
      } text ("Tab-delimited file containing old gene ID -> new gene ID representing gene-families.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.descriptionFile)
      verifyFile(config.stringsFile)
      verifyFile(config.geneMappingFile)
      fetchWithStrings(config)
    }
  }

  def fetchWithStrings(config: Config): Unit = {
    //load strings
    val strings = openFileWithIterator(config.stringsFile).toList
    println(timeStamp + "Found " + strings.size + " strings to parse")
    //load descriptions into list of id -> set(strings)
    val descriptions = loadDescriptions(config.descriptionFile)
      .map(x => (x.id, x.description.toLowerCase().split("[\\s+|\\W+]").filter(_.nonEmpty).toSet))
    println(timeStamp + "Found " + descriptions.size + " gene descriptions")
    //iterate through each string and fetch genes with keys
    val corresponding_genes = strings.foldLeft(List[(Int,String)]())((retained, string) => {
      //parse out comma separated keys
      val keys = string.toLowerCase().split("[\\s+|\\W+]").filter(_.nonEmpty)
      //retain all genes with all keys parsed and update retained list
      descriptions.filter(d => keys.forall(k => d._2(k))).foldLeft(retained)((acc, desc) => (desc._1, string) :: acc)
    })
    println(timeStamp + "Found " + corresponding_genes.size + " matching genes")
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    //load mappings, if they exist
    val mappings = {
      if(config.geneMappingFile == null) Map.empty[Int,Int]
      else
        openFileWithIterator(config.geneMappingFile).toList.map(x => {
          val c = x.split("\t")
          c(0).toInt -> c(1).toInt
        }).toMap
    }
    if(config.geneMappingFile != null) println(timeStamp + "Found " + mappings.size + " gene mappings")
    corresponding_genes.map(x =>(mappings.getOrElse(x._1, x._1), x._2)).distinct.sortBy(_._1)
      .foreach(x => pw.println(x._1 + "\t" + x._2))
    pw.close()
    println(timeStamp + "Successfully completed!")
  }

}
