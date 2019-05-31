package misc

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.GFFutils.parseMultiGFF2FingerTree
import utilities.SequenceUtils.{dna2Protein, reverseComplement}

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 16-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object ExtractGenes {

  case class Config(
                     allGenomes: File = null,
                     allGFFs: File = null,
                     outputDir: File = null,
                     featureType: String = "CDS",
                     nameTag: String = "product=",
                     idTag: String = "protein_id=",
                     prefix: String = null)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("extract-genes") {
      opt[File]('g', "genomes") required() action { (x, c) =>
        c.copy(allGenomes = x)
      } text ("File containing full path to genome assemblies, one per line.")
      opt[File]('a', "annotations") required() action { (x, c) =>
        c.copy(allGFFs = x)
      } text ("File containing full path to GFF3-formatted files, one per line.")
      opt[String]('p', "prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output files.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      note("\nOPTIONAL")
      opt[String]("feature-type") action { (x, c) =>
        c.copy(featureType = x)
      } text ("Feature types in the GFF3 file to analyse (default is 'CDS').")
      opt[String]("name-tag") action { (x, c) =>
        c.copy(nameTag = x)
      } text ("Name tag storing gene name description in the attributes column (default is 'product=')")
      opt[String]("id-tag") action { (x, c) =>
        c.copy(idTag = x)
      } text ("Name tag storing gene ID (default is 'protein_id=')")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.allGenomes)
      verifyFile(config.allGFFs)
      extractGenes(config)
    }
  }

  def extractGenes(config: Config): Unit = {
    //load all genomes
    val genomes = openFileWithIterator(config.allGenomes).toList.map(new File(_))
    //verify that they exist
    genomes.foreach(verifyFile(_))
    println(timeStamp + "Found " + genomes.size + " genome assemblies")
    //load all annotations
    val annotations = openFileWithIterator(config.allGFFs).toList.map(new File(_))
    //verify that they exist
    annotations.foreach(verifyFile(_))
    println(timeStamp + "Found " + annotations.size + " gene annotations")
    println(timeStamp + "Parsing annotations")
    //parse GFF3 files
    val (seq2Coords, id2Desc) = {
      val tmp = parseMultiGFF2FingerTree(Set(config.featureType), config.nameTag, config.idTag)(annotations, false)
      (tmp._1.mapValues(_.toList.map(x => (x._2.id, x._1, x._2.ori)).sortBy(_._1)), tmp._2)
    }
    assert(seq2Coords.values.foldLeft(0)((acc,l) => acc + l.size) == id2Desc.size)
    println(timeStamp + "--Found " + id2Desc.size + " annotations")
    //create seq output file
    val pw_seqs = new PrintWriter(config.outputDir + "/" + config.prefix + ".fna")
    val pw_proteins = new PrintWriter(config.outputDir + "/" + config.prefix + ".faa")
    println(timeStamp + "Extracting sequences")
    //iterate through genomes and output extracted genes
    genomes.foreach(genome => {
      //load iterator
      val iterator = openFileWithIterator(genome)

      /**
        * Function to
        *
        * @param seqID
        * @param seq
        */
      def outputSubSeqs(seqID: String, seq: String): Unit = {
        val coords = seq2Coords.get(seqID)
        //seq size
        val seq_size = seq.size
        if (coords.nonEmpty) {
          coords.get.foreach { case (id, (start, end), ori) => {
            //only continue if coordinates are met
            if(end <= seq_size) {
              //set gene dna sequence
              val geneseq = seq.substring(start, end).toUpperCase
              //set protein seq
              val proteinseq = dna2Protein(if(ori == '+') geneseq else reverseComplement(geneseq))
              pw_proteins.println(">" + id + "\n" + proteinseq)
              pw_seqs.println(">" + id + "\n" + geneseq)
            }
          } }
        }
      }

      /**
        * Tail-recursive method to iterate through sequence file and output specified (sub)sequences in the
        * configuration map
        *
        * @param name    Name of fasta entry as option
        * @param acc_seq Accumulating sequence as string builder
        */
      @tailrec def parseFastaEntry(name: Option[String], acc_seq: StringBuilder): Unit = {
        //iterator is empty, attempt to fetch last entry
        if (iterator.isEmpty) outputSubSeqs(name.get, acc_seq.toString)
        else {
          //get next line
          val current_line = iterator.next()
          //started new fasta entry
          if (current_line.startsWith(">")) {
            //there is an existing fasta entry, output that to disk if it was specified
            if (name.nonEmpty) outputSubSeqs(name.get, acc_seq.toString)
            //start building new fasta entry
            parseFastaEntry(Option(parseSeqName(current_line)), new StringBuilder())
          }
          //still in the same fasta entry, update sequence
          else parseFastaEntry(name, acc_seq.append(current_line))
        }
      }

      parseFastaEntry(None, StringBuilder.newBuilder)
    })
    val pw_desc = new PrintWriter(config.outputDir + "/" + config.prefix + ".gene_description.txt")
    //iterate through description and output
    id2Desc.toList.sortBy(x => (x._2.ref, (x._2.start, x._2.end)))
      //iterate through each gene and output
      .foreach(x => pw_desc.println(x._2.toString()))
    pw_desc.close()
    pw_seqs.close()
    pw_proteins.close()
    println(timeStamp + "Successfully completed!")
  }

  /**
    * Function to parse name of a sequence (fasta/fastq) entry. Only consider first string after splitting for any
    * whitespace
    *
    * @return String
    */
  private def parseSeqName: String => String = name => name.split("\\s+").head.substring(1)


}
