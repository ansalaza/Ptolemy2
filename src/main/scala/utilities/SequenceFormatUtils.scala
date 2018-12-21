package utilities

import java.io.File

import utilities.FileHandling.{getFileExtension, openFileWithIterator}
import atk.ProgressBar.progress

/**
  * Author: Alex N. Salazar
  * Created on 10-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

object SequenceFormatUtils {

  /**
    * Abstract class for Fasta and Fastq sequence entries
    */
  abstract class Entry {
    val name: String
    val seq: String
  }

  /**
    * Abstract class for Fasta and Fastq sequence iterators
    */
  abstract class SeqIterator {
    val iterator: Iterator[String]

    def hasNext(): Boolean = iterator.hasNext

    def next(): Entry
  }

  /**
    * Function to get file type for a given file
    *
    * @return String
    */
  def getFileType: File => String = file => {
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
    * Function to load generic sequence iterator
    * @return SeqIterator
    */
  def loadSeqIterator: File => SeqIterator = file => {
    //get file type
    val file_type = getFileType(file)
    //load seq iterator
    file_type match {
      case "fasta" => new FastaIterator(file)
      case "fastq" => new FastqIterator(file)
    }
  }

  /**
    * Function to construct map of sequence ID -> sequence length
    * @return Map[String, Int]
    */
  def computeSeqLengthMap: File => Map[String, Int] = reads_file => {
    //load iterator
    val iterator = loadSeqIterator(reads_file)

    /**
      * Method to iterate through seq iterator and add read name and length to list
      * @param acc Accumulating list of read name and length
      * @return List[(String,Int)]
      */
    def appendLengths(acc: List[(String, Int)]): List[(String, Int)] = {
      if (!iterator.hasNext()) acc
      else {
        //get next fasta entry
        val next_entry = iterator.next()
        //add to list
        appendLengths((next_entry.name, next_entry.seq.length) :: acc)
      }
    }
    //get read length list and covert to map data structure
    appendLengths(List()).toMap
  }

  /**
    * Method to iterate through a file of sequences (e.g reads) and assign unique IDs
    * @return
    */
  def assignSeqIDs: File => Map[String,Int] = reads_file => {
    //load iterator
    val iterator = loadSeqIterator(reads_file)

    /**
      * Method to iterate through seq iterator and assign unique INT ID to read name
      * @param acc Accumulating list of 2-tuples as (read name, unique ID)
      * @param id Accumulating ID
      * @return List[(String, Int)]
      */
    def assignIDs(acc: List[(String, Int)], id: Int): List[(String, Int)] = {
      if (!iterator.hasNext()) acc
      else {
        //get next fasta entry
        val next_entry = iterator.next()
        //add to list
        assignIDs((next_entry.name, id) :: acc, id+1)
      }
    }
    //iterate and assignd unique INT ids
    assignIDs(List(), 0).toMap
  }

  /**
    * Case class for a Fasta entry sequence
    */
  case class FastaEntry(name: String, seq: String) extends Entry

  /**
    *
    * @param file
    */
  class FastaIterator(file: File, logger: Int = 10000) extends SeqIterator {
    //load file as iterator
    val iterator = openFileWithIterator(file)
    //create an empty variable storing the current fasta entry name
    private var entry_acc = Option.empty[String]

    /**
      * Check that the iterator is not empty
      *
      * @return Boolean
      */
    //def hasNext(): Boolean = iterator.hasNext

    /**
      * Get the next Fasta entry sequence
      *
      * @return
      */
    def next(): FastaEntry = {
      progress(logger)
      /**
        * Method to fetch the next Fasta entry sequence
        *
        * @param seq_acc Accumulating sequence
        * @return FastaEntry
        */
      def buildFastaEntry(seq_acc: StringBuilder): FastaEntry = {
        /**
          * Method to create a FastaEntry object using current data accumulated
          *
          * @return FastaEntry
          */
        def createFastaEntry(): FastaEntry = {
          assert(seq_acc.nonEmpty, "Start of new Fasta entry but current entry sequence is empty " + entry_acc.get)
          //create fasta entry
          val fasta_entry = new FastaEntry(entry_acc.get.split("\\s++").head.substring(1), seq_acc.toString)
          //return
          fasta_entry
        }

        //if empty assume this is the first iteration
        if (entry_acc.isEmpty) {
          //get entry name
          val entry_name = iterator.next()
          assert(entry_name.startsWith(">"), "Expected start of new entry but found the following line: " + entry_name)
          //set current fasta entry to current entry name
          entry_acc = Option(entry_name)
          //move on to next line
          buildFastaEntry(seq_acc)
        }
        //assume building sequence
        else {
          //reached end of file, create fasta entry
          if (!hasNext()) createFastaEntry()
          else {
            //get next entry sequence
            val entry_seq = iterator.next()
            //current line is a sequence line, append to current accumulating sequence
            if (!entry_seq.startsWith(">")) buildFastaEntry(seq_acc.append(entry_seq))
            //start of new fasta entry
            else {
              val fasta_entry = createFastaEntry()
              //store next fasta entry name
              entry_acc = Option(entry_seq)
              //return
              fasta_entry
            }
          }
        }
      }
      //attempt to create the next fasta entry
      buildFastaEntry(new StringBuilder())
    }
  }

  /**
    *
    * @param name
    * @param seq
    */
  case class FastqEntry(name: String, seq: String, qual: String) extends Entry

  /**
    *
    * @param file
    */
  class FastqIterator(file: File, logger: Int = 10000) extends SeqIterator {
    //load file as iterator
    val iterator = openFileWithIterator(file)

    /**
      * Check that the iterator is not empty
      *
      * @return Boolean
      */
    //def hasNext(): Boolean = iterator.hasNext

    /**
      * Get the next Fasta entry sequence
      *
      * @return
      */
    def next(): FastqEntry = {
      progress(logger)
      /**
        * Method to iterate through each line and gather all lines for the next fastq entry
        *
        * @param acc Accumulating list in reverse order
        * @return FastqEntry
        */
      def buildFastqEntry(acc: List[String]): FastqEntry = {
        //assume each fastq entry has 4 lines, get next line and add to acc
        if (acc.size != 4) buildFastqEntry(iterator.next() :: acc)
        else {
          //get name, sequence, and corresponding quality scores
          val (name, seq, qual) = {
            //reverse list
            val tmp = acc.reverse
            //return assumed elements
            (tmp.head, tmp(1), tmp(3))
          }
          assert(name.startsWith("@"), "Unexpected start of fastq entry " + name)
          assert(seq.length == qual.length, "Unexpected quality and sequence length for fastq entry " + name)
          //return fastq entry
          new FastqEntry(name.split("\\s++").head.substring(1), seq, qual)
        }
      }

      buildFastqEntry(List())
    }
  }

}
