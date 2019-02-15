package genome_alignment

import java.io.File

import utilities.FileHandling.{verifyDirectory, verifyFile, openFileWithIterator, timeStamp, getFileName}

import utilities.Minimap2Utils.pairwiseAlignment
import utilities.SequenceFormatUtils.computeSeqLengthMap

/**
  * Author: Alex N. Salazar
  * Created on 10-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object PairwiseAlignment {

  case class Config(
                     indexDir: File = null,
                     genomesFile: File = null,
                     outputDir: File = null
                   )

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("pairwise-alignment") {
      opt[File]('g', "genomes") required() action { (x, c) =>
        c.copy(genomesFile = x)
      } text ("File containing full path to genome assemblies, one per line.")
      opt[File]('i', "index-directory") required() action { (x, c) =>
        c.copy(indexDir = x)
      } text ("Directory containing only the indexes of genomes to be alignment (using minimap2).")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyDirectory(config.indexDir)
      pairwiseGenomeAlignment(config)
    }
  }

  def pairwiseGenomeAlignment(config: Config): Unit = {
    //load genome indeces
    val all_indeces = config.indexDir.listFiles().toList.sortBy(_.getName)
    //get all genome assembly paths
    val all_assemblies = openFileWithIterator(config.genomesFile).toList.map(new File(_)).sortBy(_.getName)
    //zip together indeces and assemblies and sort by genome length
    val assembly2index = {
      //sanity check that number of files match
      assert(all_assemblies.size == all_indeces.size, "Number of assemblies does not match number of index files")
      //zip assembly to file
      val tmp = all_assemblies.zip(all_indeces)
      //sanity check
      tmp.foreach{case(assembly,index) => {
        //verify file exist
        verifyFile(assembly)
        //verify file exist
        verifyFile(index)
        //verify they correspond
        assert(getFileName(assembly) == getFileName(index), "Name of assembly file does not match index file")
      }}
      //sort by largest genome
      tmp.sortBy(x => -computeSeqLengthMap(x._1).values.sum)
    }
    println(timeStamp + "Found " + assembly2index.size + " genomes")
    println("Performing pairwise genome alignment:")
    pairwiseAlignment(assembly2index, config.outputDir)
    println("Successfully completed!")
  }

}
