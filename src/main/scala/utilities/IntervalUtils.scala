package utilities

import de.sciss.fingertree.RangedSeq
import utilities.NumericalUtils.{max, min}

import scala.annotation.tailrec

/**
  * Author: Alex N. Salazar
  * Created on 17-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object IntervalUtils {

  /**
    * Implicit way to get gene names from finger tree
    *
    * @param it
    */
  implicit class Names(it: Iterator[(_, _)]) {
    def names = it.map(_._2).mkString(", ")
  }

  /**
    * Size of an interval on in respects to a DNA/RNA sequence (+1 added).
    * @return Int
    */
  def intervalSize: ((Int,Int)) => Int = interval => (interval._2 - interval._1) + 1

  /**
    * Compute the left and right-most boundaries of a given list of intervals
    * @return (Int,Int)
    */
  def mergeIntervals: List[(Int,Int)] => (Int,Int) = intervals => (intervals.minBy(_._1)._1, intervals.maxBy(_._2)._2)

  /**
    * Method to compute the overlap of intervals. A value of at least zero means there is an overlap
    * @param x Interval one
    * @param y Interval two
    * @return Int
    */
  def computeOverlap(x: (Int,Int), y: (Int,Int)): Int = min(x._2, y._2) - max(x._1, y._1)

  /**
    * Method to compute the intersection between two intervals, if it exist. Here, intersection requires at least
    * more than 1nt
    * @param x Interval
    * @param y
    * @return
    */
  def computeIntersection(x: (Int,Int), y: (Int,Int)): Option[(Int,Int)] = {
    val start = max(x._1, y._1)
    val end = min(x._2, y._2)
    //require intersection to be more than 1 nt
    if(end - start <= 0) None else Option((start,end))
  }

  /**
    * Method to determine the alignment coverage of a gene given the alignment coordinates
    * @param gene_range
    * @param read_position
    * @return Double
    */
  def rangeCoverage(gene_range: (Int, Int), read_position: (Int,Int)): Double = {
    //overlap of the two orfs reaches the threshold
    val overlap = computeOverlap(gene_range, read_position)
    //compute size of gene
    val gene_size = intervalSize(gene_range).toDouble
    //compute coverage
    overlap/gene_size
  }

  /**
    * Method to identify all longest overlapping intervals and merge individual instances into a single interval.
    * Returns a list of sorted intervals representing merged instances of longest overlapping intervals.
    * @return List[(Int,Int)]
    */
  def longestOverlappingIntervals: List[(Int,Int)] => List[(Int,Int)] = intervals => {
    /**
      * Tail-recursive method to obtain the longest overlapping intervals.
      * @param remaining_intervals List of intervals to process, in sorted order
      * @param acc Accumulating overlapping intervals
      * @param longest List of merged longest overlapping intervals
      * @return List[(Int,Int)]
      */
    @tailrec def _longestOverlappingIntervals(remaining_intervals: List[(Int,Int)],
                                acc: List[(Int,Int)],
                                longest: List[(Int,Int)]): List[(Int,Int)] = {
      remaining_intervals match {
          //no more intervals to process
        case Nil => if(acc.isEmpty) longest else mergeIntervals(acc) :: longest
          //intervals to process
        case head :: tail => {
          //add head to acc
          if(acc.isEmpty) _longestOverlappingIntervals(tail, head :: acc, longest)
          else{
            //check whether the current interval overlaps with any of the intervals in acc
            val overlap_exist = acc.exists(x => computeOverlap(x, head) >= 0)
            //interval overlaps
            if(overlap_exist) _longestOverlappingIntervals(tail, head :: acc, longest)
            //merge curated overlaps
            else _longestOverlappingIntervals(tail, List(head), mergeIntervals(acc) :: longest)
          }
        }
      }
    }
    //sort and run
    _longestOverlappingIntervals(intervals.sortBy(identity), List(), List()).reverse
  }

}
