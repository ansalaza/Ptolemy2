package utilities

import java.io.File

import utilities.AlignmentUtils.Alignment
import utilities.FileHandling.openFileWithIterator
import utilities.IntervalUtils.{computeOverlap, intervalSize, mergeIntervals}
import utilities.NumericalUtils.{max, min}

/**
  * Author: Alex N. Salazar
  * Created on 13-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */

object PAFutils {


  /**
    * Class for individual PAF-entry lines
    * @param qname
    * @param query_length
    * @param qcoords
    * @param rname
    * @param rcoords
    * @param mapq
    */
  case class PAFentry(qname: String,
                      query_length: Int,
                      qcoords: (Int,Int),
                      ori: Char,
                      rname: String,
                      ref_length: Int,
                      rcoords: (Int,Int),
                      mapq: Int)


  /**
    * Function to convert a string to a PAFentry object
    * @return PAFentry
    */
  def toPAFentry: String => PAFentry = line => {
    //split line
    val split = line.split("\t")
    //get orientation column
    val orientation = split(4).head
    new PAFentry(split.head, split(1).toInt, (split(2).toInt, split(3).toInt),orientation, split(5), split(6).toInt,
      (split(7).toInt, split(8).toInt), split(11).toInt)
  }

  /**
    * Simple iterator to iterate through a PAF-alignment file. Using next() returns PAFentry objects
    * @param file PAFiterator
    */
  class PAFiterator(file: File) {
    //load file as iterator
    val iterator = openFileWithIterator(file)
    //standard has next method
    def hasNext(): Boolean = iterator.hasNext
    //load next paf entry
    def next(): PAFentry = toPAFentry(iterator.next())
  }

  /**
    * Method to curate all alignments per sequence given a PAF-formatted file. Returns a map containing assigned unique
    * id to a seq name along with seq length and a list of containing all curated alignments of a seq sorted by
    * read coordinates. Optionally provide a map containing previously assigned IDs to seq names: Map[Name, (Id,
    * length)]
    * @param file PAF-formatted file
    * @return (Map[Int,(String,Int)], List[(Int, List[Alignment])])
    */
  def curateAlignmentsPerSeq(file: File,
                              minMapq: Int,
                              initial_id: Int = 0,
                              initial_map: Map[String,(Int,Int)] = Map(),
                             ): (Map[Int, (String,Int)], List[(Int, List[Alignment])]) = {
    //load paf iterator
    val paf_iterator = new PAFiterator(file)

    /**
      * Function to convert a list of paf entries (along with assigned id) to a list of alignments
      * @return List[Alignment]
      */
    def pafentries2Alginments: List[PAFentry] => List[Alignment] = entries => {
      entries.map(x => new Alignment(x.qcoords, x.rname, x.ori, x.rcoords, x.mapq))
    }


    def _curateAlignmentsPerSeq(acc_alignments: List[PAFentry],
                                 id: Int,
                                 readName2IdLength: List[(String, (Int,Int))],
                                 curated_alignments: List[(Int, List[Alignment])],
                                isReadAlignments: Boolean = false
                                ): (Map[Int, (String,Int)], List[(Int, List[Alignment])]) = {
      /**
        * Function to update read ID list given an id, name, and length, if needed
        * @return List[(Int, (String,Int))]
        */
      def updateReadIdList: (String, Int) =>  List[(String, (Int,Int))] = (name, length) => {
        //check if the read has bee previously assigned
        val previous_assigned = initial_map.get(name)
        //read has been previously assigned, add existing entry to list
        if(previous_assigned.nonEmpty) readName2IdLength
        //read has never been previously assigned, add assigned id to list
        else (name, (id, length)) :: readName2IdLength
      }

      /**
        * Function to restructure provided list to map of id -> (read name, length)
        * @return Map[Int, (String,Int)]
        */
      def returnId2NameLength: List[(String, (Int,Int))] => Map[Int, (String,Int)] = read_name_list =>
        read_name_list.map(x => (x._2._1, (x._1, x._2._2))).toMap

      //reached end of alignment ifle
      if(!paf_iterator.hasNext()){
        //no more accumulated alignments to process
        if(acc_alignments.isEmpty) (returnId2NameLength(readName2IdLength), curated_alignments)
        //add final set of read alignments
        else {
          //check if there is an existing id for current read
          val existing_id = initial_map.get(acc_alignments.head.qname)
          //update list, if needed
          (returnId2NameLength(updateReadIdList(acc_alignments.head.qname, acc_alignments.head.query_length)),
            //add alignment id to list
            (if(existing_id.nonEmpty) existing_id.get._1 else id,
              pafentries2Alginments(acc_alignments).sortBy(_.qcoords)) :: curated_alignments)
        }
      }
      //still alignments to process
      else {
        //get next alignment
        val current_alignment = paf_iterator.next()
        //current alignment does not meet minimum mapping quality
        if(current_alignment.mapq < minMapq) _curateAlignmentsPerSeq(acc_alignments, id, readName2IdLength, curated_alignments)
        //add alignment to acc if it is empty (i.e. just processed a read of start of iteration)
        else if(acc_alignments.isEmpty) {
          _curateAlignmentsPerSeq(current_alignment :: acc_alignments, id, readName2IdLength, curated_alignments)
        }
          //still processing alignments from the same read
        else if(acc_alignments.head.qname == current_alignment.qname) {
          _curateAlignmentsPerSeq(current_alignment :: acc_alignments, id, readName2IdLength, curated_alignments)
        }
          //processing a new read, add accumulated alignments data structures
        else {
          //check if there is an existing id for current read
          val existing_id = initial_map.get(acc_alignments.head.qname)
          _curateAlignmentsPerSeq(
            //add current alignments to own list
            List(current_alignment),
            //update id if never seen before
            if(existing_id.nonEmpty) id else id + 1,
            //update read id list, if needed
            updateReadIdList(acc_alignments.head.qname, acc_alignments.head.query_length),
            //add curated alignment
            (if(existing_id.nonEmpty) existing_id.get._1 else id,
              pafentries2Alginments(acc_alignments).sortBy(_.qcoords)) :: curated_alignments)
        }
      }
    }
    _curateAlignmentsPerSeq(List(), initial_id, initial_map.toList, List())
  }

}
