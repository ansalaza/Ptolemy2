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
                      ori: Byte,
                      rname: String,
                      ref_length: Int,
                      rcoords: (Int,Int),
                      mapq: Int)

  /**
    * Method to classify alignment type (used for assembly).
    * @param o Overhang length
    * @param r Max overhang to mapping length ratio
    * @param a Alignment as a PAFentry object
    * @return String
    */
  def getAlignmentType(o: Int, r: Double)(a: PAFentry): String = {
    //compute the overhang of the given alignment
    val overhang = min(a.qcoords._1, a.rcoords._1) + min(a.query_length - a.qcoords._2, a.ref_length - a.rcoords._2)
    //compute the length of the alignment
    val maplen = max(intervalSize(a.qcoords), intervalSize(a.rcoords))
    //internal match
    if(overhang > min(o, (maplen * r).toInt)) "internal"
    //first reasd is contained in the second one
    else if(a.qcoords._1 <= a.rcoords._1 && (a.query_length - a.qcoords._2 <= a.ref_length - a.rcoords._2) && a.qcoords._1 < o)
      "first_contained"
    //second read is contained in the first one
    else if(a.qcoords._1 >= a.rcoords._1 && (a.query_length - a.qcoords._2 >= a.ref_length - a.rcoords._2))
      "second_contained"
    //end of the first read overlaps with the beginning of the second
    else if(a.qcoords._1 > a.rcoords._1) "first_overlap"
    //start of the first read overlaps with the end of the second
    else "second_overlap"
  }


  /**
    * Function to convert a string to a PAFentry object
    * @return PAFentry
    */
  def toPAFentry: String => PAFentry = line => {
    //split line
    val split = line.split("\t")
    //get orientation column
    val orientation = if(split(4) == "+") 0.toByte else 1.toByte
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
  /**
    * Function to find and merge overlapping alignments into a single alignment representing largest contiguous
    * stretch. Internally sorts given alignments before processing
    * @return List[PAFentry]
    */
  private def mergeAlignments: List[PAFentry] => List[PAFentry] = alignments => {
    /**
      * Function to merge a given list of overlaping alignments into a single alignment
      * @return PAFentry
      */
    def mergeOverlaps: List[PAFentry] => PAFentry = x => {
      if (x.size == 1) x.head
      else {
        //get read name
        val read_name = x.head.qname
        //get read length
        val read_length = x.head.query_length
        //get mapq
        val mapq = x.head.mapq
        //get max contiguous query coordinates
        val max_query_interval = mergeIntervals(x.map(_.qcoords))
        //get all ref sequences
        val all_refs = x.map(_.rname).distinct
        //get ref intervals
        assume(all_refs.size == 1, "Multiple reference sequences for overlapping query alignments in read " +
          read_name + ": " + all_refs)
        assume(x.map(_.ori).distinct.size == 1, "Multiple orientations for overlapping query alignments in read " +
          read_name)
        //get max contiguous ref coordinates
        val max_ref_interval = mergeIntervals(x.map(_.rcoords))
        //create merged paf entry
        new PAFentry(read_name, read_length, max_query_interval, x.head.ori, all_refs.head, x.head.ref_length,
          max_ref_interval, mapq)
      }
    }

    /**
      * Tail-recursive method to merge a given list of alignments sorted by starting position in the read (query)
      * @param remaining List of sorted alignments
      * @param acc Accumulating list of overlapping alignments
      * @param merged Final merged list of overlapping alignments
      * @return List[PAFentry]
      */
    def _mergeAlignments(remaining: List[PAFentry], acc: List[PAFentry], merged: List[PAFentry]): List[PAFentry] = {
      //no more alignments to process
      if(remaining.isEmpty){
        //no accumulated alignments to process
        if(acc.isEmpty) merged.reverse
        //merge remaining accumulated alignments
        else (mergeOverlaps(acc) :: merged).reverse
      }
      //for when acc is empty (i.e. first iteration)
      else if(acc.isEmpty) _mergeAlignments(remaining.tail, List(remaining.head), merged)
      //still alignments to process
      else {
        //alignment is on the same reference
        val is_same_ref = remaining.head.rname == acc.head.rname
        //alignment is on teh same orientation
        val is_same_ori = remaining.head.ori == acc.head.ori
        //current alignment overlaps with accumulated alignments, add to acc
        if(is_same_ref && is_same_ori && computeOverlap(remaining.head.qcoords, acc.head.qcoords) > 0)
          _mergeAlignments(remaining.tail, remaining.head :: acc, merged)
        //current alignment does not overlap with accumulated alignment, process acc and add current to new acc
        else {
          //merge overlapping alignments
          val merged_alignment = mergeOverlaps(acc)
          _mergeAlignments(remaining.tail, List(remaining.head), merged_alignment :: merged)
        }
      }
    }
    _mergeAlignments(alignments.sortBy(_.qcoords._1), List(), List()).reverse
  }
  */

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
    _curateAlignmentsPerSeq(List(), 0, initial_map.toList, List())
  }

}
