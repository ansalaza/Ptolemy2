package utilities

import java.io.{File, PrintWriter}

import utilities.FileHandling.{getFileExtension, openFileWithIterator}

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 10-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

object SequenceFormatUtils {

  /**
    * Function to get file type for a given file
    *
    * @return String
    */
  private def getFileType: File => String = file => {
    //all fasta extensions
    val fasta_types = Set("fa", "fasta", "fna")
    //all fastq extensions
    val fastq_types = Set("fastq", "fq")
    //get extension of file
    val file_extension = getFileExtension(file)
    //match file type
    file_extension match {
      //fasta type
      case x if (fasta_types.contains(x)) => "fasta"
      //fastq type
      case x if (fastq_types.contains(x)) => "fastq"
      case _ => {
        assert(false, "Unexpected file type for file " + file.getName)
        ""
      }
    }
  }

  /**
    * Internal strings for classifying fasta or fastq iles
    */
  private val fasta_type = "fasta"
  private val fastq_type = "fastq"

  /**
    * Function to parse name of a sequence (fasta/fastq) entry. Only consider first string after splitting for any
    * whitespace
    * @return String
    */
  private def parseSeqName: String => String = name => name.split("\\s+").head.substring(1)


  /**
    * Function to construct map of sequence ID -> sequence length
    * @return Map[String, Int]
    */
  def computeSeqLengthMap: File => Map[String,Int] = seq_file => {
    //load file as iterator
    val iterator = openFileWithIterator(seq_file)

    /**
      * Method to obtain a map of sequence ID -> sequence length
      * @return Map[String,Int]
      */
    def fasta2LengthMap(): Map[String, Int] = {
      /**
        * Tail-recursive method to iterate through sequence iterator and construct 2-tuple list of read ID and length
        * @param name Current fasta entry as option
        * @param length Size of current fasta entry
        * @param acc Accumulating list
        * @return 2-tuple list of read ID and length
        */
      @tailrec def _fasta2LengthMap(name: Option[String], length: Int, acc: List[(String,Int)]): List[(String,Int)] = {
        //iterator is empty, add last entry to list
        if(iterator.isEmpty) (name.get, length) :: acc
        else {
          //get next line
          val current_line = iterator.next()
          //started new fasta entry
          if(current_line.startsWith(">")) {
            _fasta2LengthMap(Option(parseSeqName(current_line)), 0,
              (if(name.isEmpty) acc else (name.get, length) :: acc))
          }
            //still in the same fasta entry, update length
          else _fasta2LengthMap(name, length + current_line.size, acc)
        }
      }
      //run method and convert output to map
      _fasta2LengthMap(None, 0, List()).toMap
    }

    /**
      * Method to obtain map of read ID -> length
      * @return Map[String, Int]
      */
    def fastq2LengthMap(): Map[String, Int] = {
      //assumes each fastq entry is in groups of 4
      iterator.grouped(4).foldLeft(List[(String,Int)]())((acc, entry) => {
        //get read id
        val name = parseSeqName(entry.head)
        //append to list
        (name, entry(1).size) :: acc
      }).toMap
    }

    //get file type
    val seq_type = getFileType(seq_file)
    //obtain length map according to file type
    seq_type match {
      case `fasta_type` => fasta2LengthMap()
      case `fastq_type` => fastq2LengthMap()
      case _ => {
        assert(false, "Unrecognized file type: " + seq_file.getName)
        Map.empty[String, Int]
      }
    }
  }


  /**
    * Method that output (sub)sequences of a given sequence file in FASTA/FASTQ format given an output file and a
    * configuration map file
    * @param config_map Configuration map containing seq ID -> list of 2-tuples as (start,end)
    * @param pw Output file as printwriter
    * @param seqs Sequence file as File object
    */
  def fetchSubSeqs(config_map: Map[String, List[(Int, Int)]], pw: PrintWriter, seqs: File): Unit = {
    /**
      * Function to output the subsequences of a sequence given the sequence name, sequence, and list of coordinates
      * @return Unit
      */
    def outputSubSeqs(name: String, seq: String): Unit = {
      //attempt to get coords, if any
      val coords = config_map.get(name)
      //this sequence exists in configuration map
      if(coords.nonEmpty) {
        //this sequence has no coordinate, output entire sequence
        if (coords.get.isEmpty) pw.println(">" + name + "\n" + seq)
        else {
          //iterate through each coordinate and output subsequence
          coords.get.foreach(coord =>
            pw.println(">" + name + "_" + coord._1 + "_" + coord._2 + "\n" + seq.substring(coord._1, coord._2))
          )
        }
      }
    }

    //load file as iterator
    val iterator = openFileWithIterator(seqs)

    /**
      * Method to iterate through FASTA file and output specified (sub)sequences in the configuration map
      */
    def fetchSubSeqFasta(): Unit = {
      /**
        * Tail-recursive method to iterate through sequence file and output specified (sub)sequences in the
        * configuration map
        * @param name Name of fasta entry as option
        * @param acc_seq Accumulating sequence as string builder
        */
      @tailrec def _fetchSubSeqFasta(name: Option[String], acc_seq: StringBuilder): Unit = {
        //iterator is empty, attempt to fetch last entry
        if(iterator.isEmpty) outputSubSeqs(name.get, acc_seq.toString)
        else {
          //get next line
          val current_line = iterator.next()
          //started new fasta entry
          if(current_line.startsWith(">")) {
            //there is an existing fasta entry, output that to disk if it was specified
            if(name.nonEmpty) outputSubSeqs(name.get, acc_seq.toString)
            //start building new fasta entry
            _fetchSubSeqFasta(Option(parseSeqName(current_line)), new StringBuilder())
          }
          //still in the same fasta entry, update sequence
          else _fetchSubSeqFasta(name, acc_seq.append(current_line))
        }
      }
      //run method and convert output to map
      _fetchSubSeqFasta(None, new StringBuilder)
    }

    /**
      * Method to iterate through FASTQ file and output specified (sub)sequences in the configuration map.
      * Note, assume FASTQ entries are grouped in 4
      */
    def fetchSubSeqFastq(): Unit = iterator.grouped(4).foreach(e => outputSubSeqs(parseSeqName(e.head), e(1)))

    //get file type
    val seq_type = getFileType(seqs)
    //obtain length map according to file type
    seq_type match {
      case `fasta_type` => fetchSubSeqFasta()
      case `fastq_type` => fetchSubSeqFastq()
      case _ => assert(false, "Unrecognized file type: " + seqs.getName)
    }
  }

}
