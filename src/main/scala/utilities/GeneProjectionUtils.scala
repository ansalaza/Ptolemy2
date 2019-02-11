package utilities

import de.sciss.fingertree.RangedSeq
import utilities.AlignmentUtils.Alignment
import utilities.GFFutils.Gene
import utilities.IntervalUtils.rangeCoverage

/**
  * Author: Alex N. Salazar
  * Created on 21-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GeneProjectionUtils {
  /**
    * Type alias for finger tree data structure
    */
  type FingerTree = RangedSeq[((Int, Int), Gene), Int]

  /** Empty figer tree data structure */
  val empty_fingertree = RangedSeq.empty[((Int, Int), Gene), Int](_._1, Ordering.Int)

  /**
    * Method to normalize position  in respects to some range [a,b] using the normalization variable function:
    * {(b-a)*(x - min_x) / (max_x - min_x))} + a
    *
    * @param a     Left-boundary of range
    * @param b     Right-boundary of range
    * @param min_x Minimum value for all values in X
    * @param max_x Maximum value for all values in X
    * @param x     Value to be normalized
    * @return Normalized value of x within range [a,b]
    */
  def normalizePosition(a: Double, b: Double, min_x: Double,
                        max_x: Double, x: Double): Int = (((b - a) * (x - min_x) / (max_x - min_x)) + a).toInt

  /**
    * Project genes onto a sequence given a list of alignments and the corresponding fingertree of the target genome.
    * It is NOT orientation aware regarding the order of genes and the orientation of genes
    *
    * @param cov_filter  Function to extract overlapping genes given coordinates from the alignment and list of
    *                    overlapping genes from fingertree
    * @param alignment   List of alignments in a target genome
    * @param local_genes Fingertree of target genome
    * @return Gene projections onto a sequence as list of gene IDs
    */
  def projectGenes(cov_filter: ((Int, Int), List[((Int, Int), Gene)]) => List[((Int, Int), Gene)])
                  (alignment: Alignment,
                   local_genes: FingerTree,
                  ): List[Gene] = {
    //get alignment coordinates on reference genome
    val ref_coords = (alignment.rcoords._1, alignment.rcoords._2 + 1)
    //get local projected genes, in order
    cov_filter(ref_coords, local_genes.filterOverlaps(ref_coords).toList).map(_._2)
    //adjust order based on alignment orientation
    //(if (alignment.isForward()) projected else projected.reverse)
    //adjust gene orientation based on alignment orientation
    //.map(x => if(alignment.isForward()) x else x.reverse())
  }

  /**
    * Function to adjust the order/orientation of a list of 2-tuples (projections and their alignment orientation
    * onto a sequence) into a final projection onto the original sequence
    *
    * @return List[Gene]
    */
  def adjustProjectionOri: List[(List[Gene], Char)] => List[Gene] = projections => {
    /**
      * Tail-recursive method to adjust orientations of a given list of 2-tuples as (list of genes, corresponding
      * alignment orientation) into a final projection onto the original sequence
      *
      * @param remaining
      * @param acc
      * @param projected
      * @return
      */
    def adjust(remaining: List[(List[Gene], Char)],
               acc: (List[Gene], Char),
               projected: List[Gene]): List[Gene] = {
      /**
        * Function to concatenate a given list of genes onto the final projection, accordingly
        *
        * @return List[Gene]
        */
      def concatenateProjections(): List[Gene] = {
        //if reverse, list is already reverse, just reverse internal gene orientations; else reverse order
        projected ::: (if (acc._2 == '-') acc._1.map(_.reverse) else acc._1.reverse)
      }
      //no more projections to process
      remaining match {
        case Nil => {
          assert(acc._1.nonEmpty)
          //concatenate final projection
          concatenateProjections()
        }
        //projections remaining
        case (head :: tail) => {
          //same orientation, update acc by prepending to list
          if (head._2 == acc._2) adjust(tail, (head._1.foldLeft(acc._1)((b, a) => a :: b), acc._2), projected)
          //different orientation concatenate accordingly
          else adjust(tail, head, concatenateProjections)
        }
      }
    }

    adjust(projections.tail, (projections.head._1.reverse, projections.head._2), List())
  }

  /**
    * Curried funtion to filter out overlapping genes retrieved from the finger tree. By definition, this means
    * comparing the first and last genes and checking the alignment coverage on whether they meet the minimum
    * alignment coverage threshold
    *
    * @param minAlignmentCov Minimum alignment coverage threshold
    * @return
    */
  def filterByCoverage(minAlignmentCov: Double
                      ): ((Int, Int), List[((Int, Int), Gene)]) => List[((Int, Int), Gene)] = (read_range, genes) => {
    //no genes to process, move one
    if (genes.isEmpty) genes
    //only one gene to process
    else if (genes.size == 1) if (rangeCoverage(genes.head._1, read_range) >= minAlignmentCov) genes else genes.tail
    //multiple genes to process
    else {
      //compute alignment coverage from the left-end
      val left_end_coverage = rangeCoverage(genes.head._1, read_range)
      //compute alignment coverage from the right-end
      val right_end_coverage = rangeCoverage(genes.last._1, read_range)
      //potentially remove left-end gene
      val tmp = if (left_end_coverage >= minAlignmentCov) genes else genes.tail
      //potentially remove right-end gene and return
      if (right_end_coverage >= minAlignmentCov) tmp else tmp.dropRight(1)
    }
  }

}
