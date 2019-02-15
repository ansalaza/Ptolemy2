package utilities

import java.io.{File, PrintWriter}

import scala.sys.process.ProcessLogger
import scala.sys.process._
import utilities.FileHandling.{getFileName, verifyFile}
import atk.ProgressBar.progress

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 10-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

object Minimap2Utils {

  /**
    * Method to index a given list of genome assemblies
    *
    * @param _parameters Additional parameters for minimap2 (default is '-k 15 -w 3')
    * @param genome      Genome assembly to index
    * @param outputdir   Output directory to store indexes
    * @return File object to indexed genome file
    */
  def indexGenome(_parameters: String, genome: File, outputdir: File): File = {
    //get additional parameters
    val parameters = _parameters.split("\\s+").toSeq

    /**
      * Function to index a genome with minimap2
      */
    def indexGenomeCommand: (File, File) => Seq[String] = (genome, output_file) =>
      Seq("minimap2") ++ parameters ++ Seq("-d", output_file.getAbsolutePath, genome.getAbsolutePath)

    //get genome name
    val name = getFileName(genome)
    //get index file
    val index_file = new File(outputdir + "/" + name + ".idx")
    println(indexGenomeCommand(genome, index_file))
    //index genome
    Process(indexGenomeCommand(genome, index_file)).!
    //verify successful index
    verifyFile(index_file, "Non-existent index file for genome: " + name)
    index_file
  }

  /**
    * Method to perform pairwise genome alignments
    *
    * @param assembly2index List of index files as File objects
    */
  def pairwiseAlignment(assembly2index: List[(File,File)], outputdir: File): Unit = {
    //create output file for storing alignments
    val pw = new PrintWriter(outputdir + "/pairwise_alignments.paf")

    //progress report
    val report = if(assembly2index.size < 10) 1 else assembly2index.size / 10

    /**
      *
      * @return
      */
    def constructAlignmentCommand: (File, List[File]) => Seq[String] = (subj, targets) => {
      Seq("minimap2", "-c") ++ Seq(subj.getAbsolutePath) ++ targets.map(_.getAbsolutePath)
    }

    /**
      *
      * @param remaining_genomes
      * @return
      */
    @tailrec def _pairwiseAlignment(remaining_genomes: List[(File,File)]): String = {
      progress(report)
      remaining_genomes match {
        //no more genomes to align
        case Nil => "Done!"
        //remaining genomes to align
        case (subj_assembly, subj_index) :: tail => {
          /*
            Capture stdout and stderr.
            Note: using mutable stringbuilder, dirty but gets job done; may be optimizedlater
              */
          var out = new StringBuilder
          var err = new StringBuilder
          val logger = ProcessLogger((o: String) => out.append(o + "\n"), (e: String) => err.append(e))
          //set command
          val command = constructAlignmentCommand(subj_assembly, tail.map(_._1))
          println("Running command")
          //run alignment command
          Process(command).!(logger)
          println("Parsing output")
          //process alignments
          out.toString.split("\n").foreach(line => if (line.nonEmpty) pw.println(line))
          //return remaining genomes
          _pairwiseAlignment(tail)
        }
      }
    }
    //perform pairwise alignments
    _pairwiseAlignment(assembly2index)
    pw.close
  }
}