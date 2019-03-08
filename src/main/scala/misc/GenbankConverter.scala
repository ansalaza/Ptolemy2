package misc

import java.io.File

import utilities.FileHandling.{verifyDirectory, verifyFile, timeStamp}
import utilities.SequenceFormatUtils.convertFromGenBank

/**
  * Author: Alex N. Salazar
  * Created on 8-3-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GenbankConverter {

  case class Config(
                     genbankFile: File = null,
                     outputDir: File = null
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("genbank-converter") {
      opt[File]('g', "genbank-file") required() action { (x, c) =>
        c.copy(genbankFile = x)
      } text ("Genome assembly in genbank format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.genbankFile)
      convertGenbank(config)
    }
  }

  def convertGenbank(config: Config): Unit = {
    convertFromGenBank(config.genbankFile, config.outputDir)
    println(timeStamp + "Successfully completed!")
  }

}
