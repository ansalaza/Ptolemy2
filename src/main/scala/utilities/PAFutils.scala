package utilities

import java.io.File

import utilities.AlignmentUtils.Alignment
import utilities.FileHandling.openFileWithIterator
import utilities.GeneProjectionUtils.normalizePosition
import utilities.IntervalUtils.{computeIntersection, computeOverlap, intervalSize, mergeIntervals}
import utilities.NumericalUtils.{abs, max, mean, min}

import scala.annotation.tailrec

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
    * Perform 1D single-linkage clustering based on a given list of breakpoints using a given max distance. Returns
    * a map containing breakpoint -> canonical breakpoint
    * @param values List of breakpoints
    * @param max_dist
    * @return Map[Int,Int]
    */
  private def breakpointClustering(values: List[Int], max_dist: Int): Map[Int, Int] = {

    /**
      * Tail-recursive method to compute 1D single linkage clustering based on pre-defined max distance. Assumes
      * input list of values is sorted
      * @param remaining List of values
      * @param acc Accumulating cluster
      * @param clusters Accumulated clusters
      * @return List[List[Int]
      */
    @tailrec def singleLinkageclustering(remaining: List[Int],
                                         acc: List[Int],
                                         clusters: List[List[Int]]): List[List[Int]] = {
      remaining match {
        //no more values to process
        case Nil => if (acc.isEmpty) clusters else acc :: clusters
        //still values to process
        case (head :: tail) => {
          //clusters with current acc
          if(acc.exists(x => abs(x - head) <= max_dist)) singleLinkageclustering(tail, head :: acc, clusters)
          //forms new cluster
          else singleLinkageclustering(tail, List(head), acc :: clusters)
        }
      }
    }
    //cluster breakpoints, create map of breakpoint -> canonical breakpoint
    singleLinkageclustering(values.sorted, List(), List()).foldLeft(Map[Int, Int]())((map, cluster) => {
      //TODO: optimize canonical breakpoint to something other than mean value?
      val mean_breakpoint = mean(cluster)
      cluster.foldLeft(map)((acc, breakpoint) => {
        //attempt to get existing value for current breakpoint
        val fetch = acc.get(breakpoint)
        //sanity check
        if(fetch.nonEmpty) assert(fetch.get == mean_breakpoint, "Breakpoint overriding existing canonical breakpoint " +
          (mean_breakpoint, cluster))
        //add canonical breakpoint
        acc + (breakpoint -> mean_breakpoint)
      })
    })
  }

  /**
    * Function to curate a given list of alignments (PAFentry) into a list of list of alignments reflecting
    * instances where an interval in the query maps to multiple references
    * @param min_overlap Minimum interval size to be considered for multi-mapping
    * @return List[List[PAFentry]
    */
  private def curateWithMultiMaps(min_overlap: Int,
                                  min_dist: Int): List[PAFentry] => List[List[PAFentry]] = entries => {

    /**
      * Function to create a list of list of alignments (PAFentry) representing multimapping instances for a given
      * list of multimapping alignments
      * @return List[List[PAFentry]
      */
    def createMultiMapRepresentations: List[PAFentry] => List[List[PAFentry]] =  multimaps => {
      //collect all breakpoints and create map of breakpoint -> canonical breakpoint
      val bmap =
        breakpointClustering(multimaps.foldLeft(List[Int]())((b, a) => a.qcoords._2 :: (a.qcoords._1 :: b)), min_dist)
      //normalize alignments to list of normalized breakpoints and their ref name
      val normalized_alignments = multimaps.map(x => ((bmap(x.qcoords._1), bmap(x.qcoords._2)), x))
      //set alignment intervals
      val alignment_intervals = bmap.values.toSet.toList.sorted.sliding(2).toList.map(x => (x.head, x(1)))
      //iterate through alignment intervals, find alignments that overlap with interval, and create new PAFentry
      // representing multimapping instance
      alignment_intervals.foldLeft(List[List[PAFentry]]())((curated, interval) => {
        //get all references that overlap with current interval and create modified PAFentries
        val overlapping_refs = normalized_alignments.foldLeft(List[PAFentry]())((acc, alignment) => {
          /**
            * Function to project a given interval coordinate onto reference coordinate
            * @return Int
            */
          def newCoord: Int => Int = coord => {
            normalizePosition(
              alignment._2.rcoords._1,
              alignment._2.rcoords._2,
              alignment._2.qcoords._1,
              alignment._2.qcoords._2, coord)
          }
          //compute intersection if it exists
          val intersection = computeIntersection(alignment._1, interval)
          //only add intersection if it exists
          if (intersection.isEmpty) acc
          else {
            //project interval coordinates onto corresponding reference reference
            val new_coords = (newCoord(interval._1), newCoord(interval._2))
            //create new paf entry of overlapping interval and add to list
            new PAFentry(alignment._2.qname, alignment._2.query_length, interval, alignment._2.ori,
              alignment._2.rname, alignment._2.ref_length, new_coords, alignment._2.mapq) :: acc
          }
        })
        //add to curated
        overlapping_refs :: curated
      }).sortBy(_.head.qcoords)
    }

    /**
      * Tail-recursive method to cluster alignments (PAF entries) of a read that overlap/multimap with one another.
      *
      * @param remaining Remaining list of alignments as PAF entries
      * @param acc       Accumulating list of overlapping alignments
      * @param multimaps List of all multimap alignment clusters
      * @return List[List[PAFentry]
      */
    @tailrec def clusterMultiMaps(remaining: List[PAFentry],
                                  acc: List[PAFentry],
                                  multimaps: List[List[PAFentry]]
                                 ): List[List[PAFentry]] = {
      remaining match {
        //no more alignment, update multimaps accordingly
        case Nil => if (acc.isEmpty) multimaps else acc :: multimaps
        //still alignments to process
        case (head :: tail) => {
          //accumulator is empty
          if (acc.isEmpty) clusterMultiMaps(tail, List(head), multimaps)
          else {
            //check if current alignment overlaps with any of the
            val overlap_exists = acc.exists(x => computeOverlap(head.qcoords, x.qcoords) >= min_overlap)
            //overlaps exists, add current alignment to acc
            if (overlap_exists) clusterMultiMaps(tail, head :: acc, multimaps)
            //no overlaps exists, update multimaps and move on
            else clusterMultiMaps(tail, List(head), acc :: multimaps)
          }
        }
      }
    }
    //only one alignment, multimaps cannot exists
    if (entries.size == 1) List(entries)
    else {
      //cluster alignments for multimapping and sort by left-most coordinate within each cluster
      val clustered = clusterMultiMaps(entries, List(), List()).sortBy(_.map(_.qcoords._1).min)
      //iterate through each cluster and create multimap representations
      clustered.foldLeft(List[List[PAFentry]]())((acc, alignments) => {
        //no multimapping for this interval, add as is
        if(alignments.size == 1) alignments :: acc
        //create multimapping representation, add each interval(s) to list
        else createMultiMapRepresentations(alignments).foldLeft(acc)((a,b) => b :: a)
      }).reverse
    }
  }

  /**
    * Method to curate all alignments per sequence given a PAF-formatted file. Returns a 2-tuple (map, list). The map
    * keys are the assigned unique ID and values are 2-tuples: (seq name, seq length). The list elements are
    * 2-tuples: (ID, list of list of alignments sorted by query coordinates). The nested list structure is only
    * relevant when specifying minimum overlap length (minMultiMap) > 0, in which the alignments will be curated to
    * reflect intervals in the query that align to multiple references (uses minDist for breakpoint clustering).
    * Alignments are taken as is.
    *
    * @param file PAF alignment file.
    * @param minMapq Minimum mapping quality. Only alignments of at least this MAPQ will be considered.
    * @param minMultiMap Min multimapping (overlap) length (see description).
    * @param minDist Min distance for breakpoint clustering (see description).
    * @param initial_id Initial integer to be use to assigned unique IDs, default is 0.
    * @param initial_map Initial map where keys are original seq names and values are 2-tuples: (ID, length)
    * @return 2-tuple: (map where keys are IDs and values are 2-tuples: (name, length), list of 2-tuples: (ID list
    *         of list of alignments)
    */
  def curateAlignmentsPerSeq(file: File,
                             minMapq: Int,
                             minMultiMap: Int,
                             minDist: Int,
                             initial_id: Int = 0,
                             initial_map: Map[String, (Int, Int)] = Map(),
                            ): (Map[Int, (String, Int)], List[(Int, List[List[Alignment]])]) = {
    //load paf iterator
    val paf_iterator = new PAFiterator(file)

    /**
      * Function to convert a list of paf entries (along with assigned id) to a list of alignments sorted by query
      * coordinates
      *
      * @return List[List[Alignment]
      */
    def pafentries2Alginments: List[PAFentry] => List[List[Alignment]] = entries => {
      //use alignments as is
      if(minMultiMap < 1) {
        entries.map(x => List(new Alignment(x.qcoords, x.rname, x.ori,x.rcoords, x.mapq))).sortBy(_.head.qcoords)
      }
      //curate alignments with min multimap (overlap) length and min breakpoint distance for clustering
      else {
        curateWithMultiMaps(minMultiMap, minDist)(entries.sortBy(_.qcoords))
          .map(alignments => alignments.map(x => new Alignment(x.qcoords, x.rname, x.ori,x.rcoords, x.mapq)))
          .sortBy(_.head.qcoords)
      }
    }


    def _curateAlignmentsPerSeq(acc_alignments: List[PAFentry],
                                id: Int,
                                seqName2IdLength: List[(String, (Int, Int))],
                                curated_alignments: List[(Int, List[List[Alignment]])],
                                isReadAlignments: Boolean = false
                               ): (Map[Int, (String, Int)], List[(Int, List[List[Alignment]])]) = {
      /**
        * Function to update read ID list given an id, name, and length, if needed
        *
        * @return List[(Int, (String,Int))]
        */
      def updateReadIdList: (String, Int) => List[(String, (Int, Int))] = (name, length) => {
        //check if the read has bee previously assigned
        val previous_assigned = initial_map.get(name)
        //read has been previously assigned, add existing entry to list
        if (previous_assigned.nonEmpty) seqName2IdLength
        //read has never been previously assigned, add assigned id to list
        else (name, (id, length)) :: seqName2IdLength
      }

      /**
        * Function to restructure provided list to map of id -> (read name, length)
        *
        * @return Map[Int, (String,Int)]
        */
      def returnId2NameLength: List[(String, (Int, Int))] => Map[Int, (String, Int)] = read_name_list =>
        read_name_list.map(x => (x._2._1, (x._1, x._2._2))).toMap

      //reached end of alignment ifle
      if (!paf_iterator.hasNext()) {
        //no more accumulated alignments to process
        if (acc_alignments.isEmpty) (returnId2NameLength(seqName2IdLength), curated_alignments)
        //add final set of read alignments
        else {
          //check if there is an existing id for current read
          val existing_id = initial_map.get(acc_alignments.head.qname)
          //update list, if needed
          (returnId2NameLength(updateReadIdList(acc_alignments.head.qname, acc_alignments.head.query_length)),
            //add alignment id to list
            (if (existing_id.nonEmpty) existing_id.get._1 else id,
              pafentries2Alginments(acc_alignments)) :: curated_alignments)
        }
      }
      //still alignments to process
      else {
        //get next alignment
        val current_alignment = paf_iterator.next()
        //current alignment does not meet minimum mapping quality
        if (current_alignment.mapq < minMapq) _curateAlignmentsPerSeq(acc_alignments, id, seqName2IdLength, curated_alignments)
        //add alignment to acc if it is empty (i.e. just processed a read of start of iteration)
        else if (acc_alignments.isEmpty) {
          _curateAlignmentsPerSeq(current_alignment :: acc_alignments, id, seqName2IdLength, curated_alignments)
        }
        //still processing alignments from the same read
        else if (acc_alignments.head.qname == current_alignment.qname) {
          _curateAlignmentsPerSeq(current_alignment :: acc_alignments, id, seqName2IdLength, curated_alignments)
        }
        //processing a new read, add accumulated alignments data structures
        else {
          //check if there is an existing id for current read
          val existing_id = initial_map.get(acc_alignments.head.qname)
          _curateAlignmentsPerSeq(
            //add current alignments to own list
            List(current_alignment),
            //update id if never seen before
            if (existing_id.nonEmpty) id else id + 1,
            //update read id list, if needed
            updateReadIdList(acc_alignments.head.qname, acc_alignments.head.query_length),
            //add curated alignment
            (if (existing_id.nonEmpty) existing_id.get._1 else id,
              pafentries2Alginments(acc_alignments)) :: curated_alignments)
        }
      }
    }
    _curateAlignmentsPerSeq(List(), initial_id, initial_map.toList, List())
  }

}
