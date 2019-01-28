package metrics

import java.io.{File, PrintWriter}

import utilities.AlignmentUtils.{getBreakPoints, getUnalingedRegions, reportMultiMapping}
import utilities.PAFutils.curateAlignmentsPerSeq
import utilities.FileHandling.{getFileName, timeStamp, verifyDirectory, verifyFile}
import utilities.IntervalUtils.intervalSize
import utilities.SequenceFormatUtils.computeSeqLengthMap
import utilities.NumericalUtils.roundUp

import scala.collection.immutable.HashSet

/**
  * Author: Alex N. Salazar
  * Created on 9-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

object AlignmentMetrics {

  case class Config(
                     pafFile: File = null,
                     readFile: File = null,
                     minClipping: Int = 300,
                     minMapq: Int = 10,
                     outputDir: File = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("alignment-metrics") {
      opt[File]('p', "paf-file") required() action { (x, c) =>
        c.copy(pafFile = x)
      } text ("Alignment in PAF-format.")
      opt[File]('r', "read-file") required() action { (x, c) =>
        c.copy(readFile = x)
      } text ("Reads in FASTA or FASTQ-format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      note("OPTIONAL")
      opt[Int]("min-mapq") action { (x, c) =>
        c.copy(minMapq = x)
      } text ("Minimum mapping quality for an alignment to be processed(minimum is 10).")
      opt[Int]("min-clipping") action { (x, c) =>
        c.copy(minClipping = x)
      } text ("Minimum clipping length (default is 300).")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.pafFile)
      computeAlignmentMetrics(config)
    }
  }

  def computeAlignmentMetrics(config: Config): Unit = {
    //val readname2id_map = readname2id_map
    println(timeStamp + "Obtaining alignment metrics")
    //curate alignments
    val (readid2name, curated_alignments) = {
      println(timeStamp + "Creating read length map")
      //construct read id to (name, length) map
      val (initial_readid2namelength, last_id) = {
        //construct read length map
        computeSeqLengthMap(config.readFile).toList
          //iterate and assign unique id
          .foldLeft((List[(String,(Int,Int))](), 0)) { case ((acc_list, id), (name, length)) => {
          ((name, (id, length)) :: acc_list, id+1)
          }}
      }
      curateAlignmentsPerSeq(config.pafFile, config.minMapq, last_id, initial_readid2namelength.toMap)
    }
    println(timeStamp + "Found " + curated_alignments.size + " aligned reads (" +
      ((curated_alignments.size.toDouble / readid2name.size) * 100) + "% mapped)")
    //create output file for alignment summary
    val pw_align = new PrintWriter(config.outputDir + "/" + getFileName(config.pafFile) + ".alignment_sizes.txt")
    pw_align.println("ReadName\tRef\tStart\tEnd\tSize\tMAPQ\tAlignCov")
    //create output file for breakpoints summary
    val pw_brk = new PrintWriter(config.outputDir + "/" + getFileName(config.pafFile) + ".breakpoints.txt")
    pw_brk.println("ReadName\tRef\tPosition")
    //create output file for multimapping summary
    val pw_multi = new PrintWriter(config.outputDir + "/" + getFileName(config.pafFile) + ".multimapping.txt")
    pw_multi.println("ReadName\tOstart\tOend\tRstart\tRend\tOsize")
    //create output file for multimapping summary
    val pw_unmapped = new PrintWriter(config.outputDir + "/" + getFileName(config.pafFile) + ".unmapped_reads.txt")
    pw_unmapped.println("ReadName\tSize")
    //create output file for multimapping summary
    val pw_unaligned = new PrintWriter(config.outputDir + "/" + getFileName(config.pafFile) + ".unaligned.txt")
    pw_unaligned.println("ReadName\tStart\tEnd\tSize")
    //iterate through curated alignments and output metrics
    val mapped_reads = curated_alignments.foldLeft(HashSet[Int]()){case(observed,(read_id, alignments)) => {
      //get read name
      val (read_name, read_length) = readid2name(read_id)
      //check if there are multiple alignments
      val is_multi = alignments.size > 1
      //iterate through each read and output alignment
      alignments.foreach(alignment => {
        //get size of alignment, in the reference
        val ref_align_size = intervalSize(alignment.rcoords)
        //get size of alignment, in the read
        val read_align_size = intervalSize(alignment.qcoords)
        //output metric
        pw_align.println(read_name + "\t" + alignment.ref + "\t" + alignment.rcoords._1 + "\t" + alignment.rcoords._2 +
          "\t" + ref_align_size + "\t" + alignment.mapq + "\t" + roundUp(read_align_size.toDouble / read_length, 2))
      })
      //get breakpoints
      val breakpoints = getBreakPoints(alignments, read_length, config.minClipping)
      //iterate through each breakpoint and output
      breakpoints.foreach { case (ref, pos) => pw_brk.println(read_name + "\t" + ref + "\t" + pos) }
      //get multimapping
      val multimapping = reportMultiMapping(read_name, alignments)
      //output only if multimapping exist
      if (multimapping.nonEmpty) multimapping.foreach(pw_multi.println)
      //get unaligned regions
      val unaligned = getUnalingedRegions(alignments, read_length)
      //report unaligned region, if it exist
      if(unaligned.nonEmpty) {
        unaligned.foreach(region =>
          pw_unaligned.println(read_name + "\t" + region._1 + "\t" + region._2 + "\t" + intervalSize(region)))
      }
      observed + (read_id)
    }}
    //iterate through each read and report unmapped reads
    readid2name.toList.foreach{case(id, (name, length)) => if(!mapped_reads(id)) pw_unmapped.println(name + "\t" + length)}

    pw_brk.close
    pw_align.close
    pw_multi.close
    pw_unmapped.close
    pw_unaligned.close
    println(timeStamp + "Successfully completed!")
  }

}
