package utilities

import java.io.{File, PrintWriter}

import utilities.FileHandling.{getFileExtension, getFileName, openFileWithIterator, timeStamp}

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
  def parseSeqName: String => String = name => name.split("\\s+").head.substring(1)


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
  def fetchSubSeqs(config_map: Map[String, List[(Int, Int)]],
                   pw: PrintWriter,
                   seqs: File): Unit = {
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
          //get total sequence length
          val length = seq.size
          //iterate through each coordinate and output subsequence
          coords.get.foreach(coord =>{
            //linear extraction
            if(coord._2 <= length)
              pw.println(">" + name + "_" + coord._1 + "_" + coord._2 + "\n" + seq.substring(coord._1, coord._2))
            //circular extraction
            else {
              //compute coordinates overlapping with start (circular genome
              val start_offset = coord._2 - length
              //output dna sequence considering circular layout
              pw.println(">" + name + "_" + coord._1 + "_" + (coord._2 + start_offset))
              pw.print(seq.substring(coord._1, coord._2))
              pw.println(seq.substring(0, start_offset))
            }
          }
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

  /**
    *
    * @param file
    * @param outputdir
    */
  def convertFromGenBank(file: File, outputdir: File): Unit = {

    /**
      *
      * @return
      */
    def isDNA: Char => Boolean = nt => List('A', 'C', 'T', 'G', 'N').exists(_ == nt.toUpper)

    /**
      * Parse coordinates and orientation of a gene's coordinate line. Returns a 2-tuple representing a list of
      * coordinates and orientation. Cicular-aware, hence, list of coordinates.
      * @return 2-tuple: (List[Int], Char)
      */
    def parseCoordsLine: String => (List[Int], Char) = line => {
      //get coordinates by removing tags like 'join' and 'complement', join all coords into a list of ints
      val coords =
        line.filter(x => x.isDigit || x == '.' | x == ',').replaceAll("..", ",").split(",").toList.map(_.toInt)
      //assert either two coordinates or four
      assert(coords.size == 2 || coords.size == 4, "Unexpected number of coordinates: " + line)
      //forward orientation, extract coords
      if(!line.contains("complement")) (coords, '+')
      //entire gene is in reverse orientation
      else if(line.startsWith("complement")) (coords, '-')
      //two joined-segments, both which are in reverse orientation
      else {
        //pseudo regex line
        val reg = "complement".r
        //assert there are two complement segments
        assert(reg.findAllMatchIn(line).size == 2, "Expected two complement segments: " + line)
        //assert four coordinates
        assert(coords.size == 4, "Expected four coordinates due to two complement segments: " + line)
        //get first and second set of coordinates
        val (first, second) = (coords.take(2), coords.takeRight(2))
        //assert first coordinates are all larger than the second
        assert(second.head == 1 && first.forall(x => second.forall(_ < x)),
          "Expected second segment to begin at start of the genome: " + line)
        (coords, '-')
      }
    }

    //set fasta output
    val pw_fasta = new PrintWriter(outputdir + "/" + getFileName(file) + ".fasta")
    val pw_proteins = new PrintWriter(outputdir + "/" + getFileName(file) + ".proteins.fasta")
    //set gff output
    val pw_gff = new PrintWriter(outputdir + "/" + getFileName(file) + ".gff")

    //mutable variables for seqname and seq length
    var seqname = ""
    var length = 0
    var translation = StringBuilder.newBuilder
    var isCDS = false

    //set iterator
    val iterator = openFileWithIterator(file)

    /**
      * Function to output CDS (gene) to GFF format given the coord line, the product line, and the corresponding ID
      *
      * @return Unit
      */
    def outputGene: (List[String], Int) => Unit = (lines, id) => {
      //sanity check
      assert(lines.size == 3, "Expected 3 CDS entry lines: " + lines)
      //get coordinates and orientation
      val (coords, ori) = parseCoordsLine(lines.last)
      //set start end
      val (start,end) = (coords.head, if(coords.size == 2) coords(1) else coords(1) + coords(3))
      //get gene product description
      val description = lines.head.replaceAll("/product=", "").drop(1).dropRight(1)
      //get protein sequence
      val protein_seq = lines(1).replace("/translation=", "").drop(1).dropRight(1)
      pw_gff.println(List(seqname, "", "CDS", start, end, "", ori, "",
        "protein_id=" + (seqname + "_" + id) + ";product=" + description).mkString("\t"))
      pw_proteins.println(">" + seqname + "_" + id)
      pw_proteins.println(protein_seq)
    }

    /**
      *
      * @param acc
      * @param id
      * @return
      */
    def parseGenBank(acc: List[String], id: Int): String = {
      //end of file
      if (iterator.isEmpty) {
        //no more sequences to process
        if (acc.isEmpty) "done"
        //one more sequence to process
        else {
          outputGene(acc, id)
          "done"
        }
      }
      //more lines to process
      else {
        //get next line
        val entry = iterator.next().split("\\s+").filter(_.nonEmpty)
        //match line type accordingly
        entry.head match {
          //start sequence information
          case "LOCUS" => {
            //set sequence name
            seqname = entry(1)
            //se sequenc length
            length = entry(2).toInt
            //output genome fast info
            pw_fasta.println(">" + seqname)
            println(timeStamp + "Processing " + seqname + " of length " + length)
            parseGenBank(acc, id)
          }
          case "CDS" => {
            isCDS = true
            //no existing CDS, add as new entry
            if (acc.isEmpty) parseGenBank(entry(1) :: acc, id)
            //start of new CDS, end of old CDS
            else {
              //output old CDS
              outputGene(acc, id)
              //move on with new CDS
              parseGenBank(List(entry(1)), id + 1)
            }
          }
          case _ => {
            //genome fasta line
            if (entry.head.forall(_.isDigit) && entry.drop(1).forall(_.forall(isDNA(_)))) {
              pw_fasta.println(entry.drop(1).mkString(""))
              parseGenBank(acc, id)
            }
            //start of translation info for CDS
            else if (isCDS && entry.head.startsWith("/translation")) {
              //set start of translation flag
              assert(translation.isEmpty, "Expected translation sequence to be empty: " + translation.mkString)
              translation.append(entry.head)
              parseGenBank(acc, id)
            }
            //start of product info for CDS and end of translation
            else if (isCDS && entry.head.startsWith("/product")){
              assert(translation.nonEmpty, "Expected translation sequence to exist: " + translation.mkString)
              //set protein seq
              val protein_seq = translation.mkString
              //set translation
              translation = StringBuilder.newBuilder
              //set end of CDS
              isCDS = false
              parseGenBank(entry.mkString(" ") :: (protein_seq :: acc), id)
            }
            else if(isCDS && translation.nonEmpty) {
              translation.append(entry.head)
              parseGenBank(acc, id)
            }
            //misc info
            else parseGenBank(acc, id)
          }

        }
      }
    }

    parseGenBank(List(), 0)
    pw_gff.close
    pw_fasta.close()
    pw_proteins.close

  }

}
