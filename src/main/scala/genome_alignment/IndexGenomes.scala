package genome_alignment

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}

import utilities.Minimap2Utils.indexGenome

import atk.ProgressBar.progress

/**
  * Author: Alex N. Salazar
  * Created on 10-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object IndexGenomes {


  case class Config(
                     genomesFile: File = null,
                     outputDir: File = null,
                     parameters: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("index-genomes") {
      opt[File]('g', "genomes") required() action { (x, c) =>
        c.copy(genomesFile = x)
      } text ("File containing full path to genome assemblies, one per line.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      note("OPTIONAL\n")
      opt[String]("parameters") action { (x, c) =>
        c.copy(parameters = x)
      } text ("Minimap2 alignment parameters (default is '-k 15 -w 3').")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.genomesFile)
      pairwiseGenomeAlignment(config)
    }
  }

  def pairwiseGenomeAlignment(config: Config): Unit = {
    //set parameters
    val parameters = if(config.parameters == null) "-k 15 -w 3" else config.parameters
    //load genome assemblies
    val all_genomes = openFileWithIterator(config.genomesFile).toList.map(new File(_))
    //verify genome assemblies exists
    all_genomes.foreach(verifyFile(_))
    println(timeStamp + "Found " + all_genomes.size + " genomes")
    println(timeStamp + "Indexing genomes:")
    //iterate through each genome
    all_genomes.foreach(genome => {
      progress(5)
      //index genome
      indexGenome(parameters, genome, config.outputDir)
    })
    println("Successfully completed!")
  }


}
