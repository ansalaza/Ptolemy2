package utilities

import utilities.IntervalUtils.{computeIntersection, computeOverlap, intervalSize, longestOverlappingIntervals}
import utilities.NumericalUtils.max

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 17-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object AlignmentUtils {

  /**
    * Case class for storing an alignment
    *
    * @param qcoords
    * @param rcoords
    * @param mapq
    */
  case class Alignment(qcoords: (Int, Int),
                       ref: String,
                       ori: Char,
                       rcoords: (Int, Int),
                       mapq: Int) {
    /**
      * Alignment is in the forward-strand on the reference
      *
      * @return Boolean
      */
    def isForward(): Boolean = ori == '+'

    /**
      * Function to obtain the ref coordinates of an alignment considering the orientation of the alignment
      *
      * @return (Int,Int)
      */
    def getTrueRefPositions(): (Int, Int) = if (isForward()) rcoords else rcoords.swap
  }

  /**
    * Method to obtain breakpoint coordinates in the respective reference sequence in a given list of alignments. Here,
    * breakpoints are the locations of discontiguous alignments (i.e. end/start of two alignments next to each other
    * on the read but in unique positions in some ref seq). The minimum clipping parameter is used to determine
    * whether there is a breakpoint on the left/right-most boundaries in respects to the read (i.e. if the left and
    * right ends of a read has sequence that is not aligned and meets the minimum clipping parameter, these
    * boundaries are considered breakpoints).
    *
    * @param alignments   List of alignments (assumed to be sorted by read coordinates)
    * @param total_length Total length of read
    * @param min_clip     Minimum clipping size
    * @return
    */
  def getBreakPoints(alignments: List[Alignment], total_length: Int, min_clip: Int): List[(String, Int)] = {
    //add breakpoints from the left/right ends, if they exist
    val end_breakpoints = {
      //get left and right most boundaries
      val (min_start_alignment, max_end_alignment) = (alignments.minBy(_.qcoords._1), alignments.maxBy(_.qcoords._2))
      //compute difference at the boundaries
      val (left_diff, right_diff) = (min_start_alignment.qcoords._1, (total_length - max_end_alignment.qcoords._2) + 1)
      //add left first
      val tmp = {
        //left breakpoint is not large enough
        if (left_diff < min_clip) List()
        //left breakpoint is large enough
        else List((min_start_alignment.ref, (min_start_alignment.getTrueRefPositions())._1))
      }
      //right breakpoint is not large enough
      if (right_diff < min_clip) tmp
      //right breakpoint is large enough
      else (max_end_alignment.ref, (max_end_alignment.getTrueRefPositions())._2) :: tmp
    }
    //only one alignment so return only the breakpoints from the ends
    if (alignments.size == 1) end_breakpoints
    //multiple alignments
    else {
      //add breakpoints in between alignments; assumes alignments are already sorted
      alignments.sliding(2).foldLeft(end_breakpoints)((break_points, pair) => {
        //for easy acces
        val (alignment1, alignment2) = (pair.head, pair(1))
        //add breakpoints
        (alignment2.ref, (alignment2.getTrueRefPositions())._1) ::
          ((alignment1.ref, (alignment1.getTrueRefPositions)._2) :: break_points)
      }).sortBy(_._2)
    }
  }

  /**
    * Method to report multimapping within alignments of the same read.
    *
    * @param name       Name of read
    * @param alignments All alignments of that read
    * @return List[String]
    */
  def reportMultiMapping(name: String, alignments: List[Alignment]): List[String] = {
    /**
      * Function to create a string of a multimapped instance
      *
      * @return String
      */
    def buildReport(intersection: (Int, Int), read_coords: (Int, Int), ref_coords: (Int, Int), ref: String,
                    forward: Boolean, size: Int): String = {
      /**
        * Get true interval of based on alignment orientation given the size of the overlap
        *
        * @return (Int,Int)
        */
      def getTrueInterval: Boolean => (Int, Int) = isForward => {
        if (isForward && forward) (ref_coords._1, ref_coords._1 + size) else (ref_coords._2 - size, ref_coords._2)
      }

      //get approximate intersection coordinates on ref
      val intersect_ref = {
        //left dove tail
        if (intersection._1 == read_coords._1 && intersection._2 != read_coords._2) {
          if (forward) (ref_coords._1, ref_coords._1 + size) else (ref_coords._2 - size, ref_coords._2)
        }
        //right dove
        else if (intersection._1 != read_coords._1 && intersection._2 == read_coords._2) {
          if (forward) (ref_coords._2 - size, ref_coords._2) else (ref_coords._1, ref_coords._1 + size)
        }
        //contained
        else if (intersection._1 == read_coords._1 && intersection._2 == read_coords._2) ref_coords
        //contains
        else (ref_coords._1 + intersection._1 - 1, ref_coords._2 - (ref_coords._2 - intersection._2))
      }
      name + "\t" + intersection._1 + "\t" + intersection._2 + "\t" + ref + "\t" + intersect_ref._1 + "\t" +
        intersect_ref._2 + "\t" + size
    }
    //iterate through each alignment as subj
    alignments.foldLeft(List[String]())((overlaps, subj) => {
      //iterate through each alignmetn as target
      alignments.foldLeft(overlaps)((local_overlaps, target) => {
        //subj and target and are not the same and overlap in the read
        if (subj != target && computeOverlap(subj.qcoords, target.qcoords) > 0) {
          //obtain intersection
          val intersection = computeIntersection(subj.qcoords, target.qcoords)
          assert(intersection.nonEmpty, "Expected overlapping alignments for read " + name)
          //build string report and add
          buildReport(intersection.get, subj.qcoords, subj.rcoords, subj.ref, subj.isForward(),
            intervalSize(intersection.get)) :: local_overlaps
        }
        //subj are the same alignment or there is no overlap
        else local_overlaps
      })
    })
  }


  /**
    * Method to obtain unaligned regions given a list of alignments and a total sequence length
    *
    * @param alignments   List of alignments of some sequence
    * @param total_length Total lengt of sequence
    * @return List[(Int,Int)]
    */
  def getUnalingedRegions(alignments: List[Alignment], total_length: Int): List[(Int, Int)] = {
    //construct longest overlapping intervals
    val longest_overlapping_intervals = longestOverlappingIntervals(alignments.map(_.qcoords))
    //determine whether left end of read is unaligned
    val left_unaligned = {
      //get unaligned size
      val size = longest_overlapping_intervals.head._1
      //if it meets mininum threshold report
      if (size == 0) None else Option((1, size))
    }
    //determine whether right-end of read is unaligned
    val right_unaligned = {
      //get unalinged size
      val position = longest_overlapping_intervals.last._2 + 1
      val size_diff = total_length - longest_overlapping_intervals.last._2
      if (size_diff == 0) None else Option((position, total_length))
    }

    /**
      * Function to add the unaligned regions to the front and end of a list
      *
      * @return List[(Int,Int)]
      */
    def addEndUnaligned: List[(Int, Int)] => List[(Int, Int)] = intervals => {
      val tmp = (if (left_unaligned.isEmpty) intervals else left_unaligned.get :: intervals)
      if (right_unaligned.isEmpty) tmp else tmp.:+(right_unaligned.get)
    }

    //only one alignment, add end clips if they exist
    if (longest_overlapping_intervals.size == 1) addEndUnaligned(List())
    else {
      //iterate through each interval and compute unaligned regions
      val unaligned_intervals = {
        longest_overlapping_intervals.sliding(2).foldLeft(List[(Int, Int)]())((unaligned, intervals) => {
          //get left and right interval for easy access
          val (left_interval, right_interval) = (intervals.head, intervals(1))
          assert(right_interval._1 > left_interval._2, "Expected non-overlapping intervals in unalinged regions: " +
            alignments)
          //size difference between interval
          val size_diff = (right_interval._1 - left_interval._2) - 1
          //move on if the size difference does not exist else add and continue
          if (size_diff == 0) unaligned else (left_interval._2 + 1, right_interval._1 - 1) :: unaligned
        }).reverse
      }
      addEndUnaligned(unaligned_intervals)
    }
  }

  /**
    * Function to fetch multimapping alignments given a list of alignments for a read
    *
    * @return List[List[Alignment]
    **/
  def getMultiMapping(alignments: List[Alignment], min_overlap: Double): List[List[Alignment]] = {
    /**
      * Tail-recursive method to identify multimapping alignments
      *
      * @param remaining_alignments
      * @param acc
      * @param multimapped
      * @return
      */
    @tailrec def _getMultiMapping(
                                   remaining_alignments: List[Alignment],
                                   acc: List[Alignment],
                                   multimapped: List[List[Alignment]]): List[List[Alignment]] = {
      //check whether there are still alignments to process
      remaining_alignments match {
        //no alignments to process
        case Nil => if (acc.size < 2) multimapped else acc :: multimapped
        //alignments to process
        case current_alignment :: tail => {
          //acc is empty (i.e. first iteration) add and continue
          if (acc.isEmpty) _getMultiMapping(tail, List(current_alignment), multimapped)
          else {
            //compute size of overlap between current alignment and head alignment in acc
            val overlap = computeOverlap(current_alignment.qcoords, acc.head.qcoords)
            //compute minimum alignment size
            val min_align_size = max(intervalSize(current_alignment.qcoords), intervalSize(acc.head.qcoords))
            // overlap is sufficient to be considered multi-mapped
            if (overlap >= 0 && (overlap / min_align_size.toDouble) >= min_overlap)
              _getMultiMapping(tail, current_alignment :: acc, multimapped)
            // there is no overlap and there is multimapping in the current acc
            else if (acc.size > 1) _getMultiMapping(tail, List(current_alignment), acc :: multimapped)
            // there is no overlap and there is NO multimapping in the current acc
            else _getMultiMapping(tail, List(current_alignment), multimapped)
          }
        }
      }
    }
    //only call internal method if there are at least two different alignments
    if (alignments.size < 2) List() else _getMultiMapping(alignments.sortBy(_.qcoords), List(), List())
  }

  /**
    * Method to create indirect pairwise alignments of a given multimapped region. Assumes that it is the same query
    * sequence and coordinates aligned to multiple refs. Indirectly, the coordinates of these refs and combined
    * together pairwise into an artificial alignment.
    *
    * @param multimapped_alignments List of multimapped alignments
    * @return List of 2-tuples (ID of sequence, list of indirect alignments for current sequence)
    */
  def getCanonicalMultimappedRegions(multimapped_alignments: List[Alignment]): Map[String, List[(Int,Int)]] = {

    def createCanonicalCoord: List[Alignment] => (String, (Int,Int)) = alignments => {
      (alignments.head.ref, longestOverlappingIntervals(alignments.map(_.rcoords)).head)
    }

    /**
      * Tail-recursive method to create indirect pairwise alignments between refs of the same query sequence
      *
      * @param remaining Remaining list of alignments
      * @return same as parent method
      */
    @tailrec def canonicalRegions(remaining: List[Alignment],
                                  acc: List[Alignment],
                                  canonicals: List[(String,(Int,Int))]
                                 ): List[(String,(Int,Int))] = {
      remaining match {
        case Nil => if(acc.isEmpty) canonicals else createCanonicalCoord(acc) :: canonicals
        case (head :: tail) => {
          if(acc.isEmpty) canonicalRegions(tail, List(head), canonicals)
          else {
            //overlaps with current acc
            if(head.ref == acc.head.ref && acc.exists(a => computeOverlap(head.rcoords, a.rcoords) > 0))
              canonicalRegions(tail, head :: acc, canonicals)
            else canonicalRegions(tail, List(head), createCanonicalCoord(acc) :: canonicals)
          }
        }
      }
    }
    println("HERE")
    multimapped_alignments.sortBy(x => (x.ref, x.rcoords)).foreach(println)
    //create all pairwise instances, group by seq ID, create list of 2-tuples
    canonicalRegions(multimapped_alignments.sortBy(x => (x.ref, x.rcoords)), List(), List())
      .groupBy(_._1).mapValues(_.map(_._2))
  }

}
