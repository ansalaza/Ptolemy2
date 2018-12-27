package utilities

import de.sciss.fingertree.RangedSeq
import utilities.AlignmentUtils.Alignment
import utilities.GeneGraphUtils.PathEntry
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
  type FingerTree = RangedSeq[((Int, Int), Int), Int]

  /** Empty figer tree data structure */
  val empty_fingertree = RangedSeq.empty[((Int, Int), Int), Int](_._1, Ordering.Int)

  /**
    * Project genes onto a sequence given a list of alignments and the corresponding fingertree of the target genome.
    * It is orientation aware.
    *
    * @param cov_filter  Function to extract overlapping genes given coordinates from the alignment and list of
    *                    overlapping genes from fingertree
    * @param alignment   List of alignments in a target genome
    * @param local_genes Fingertree of target genome
    * @return Gene projections onto a sequence as list of gene IDs
    */
  def projectGenes(cov_filter: ((Int, Int), List[((Int, Int), Int)]) => List[((Int, Int), Int)])
                  (alignment: Alignment,
                   local_genes: FingerTree,
                  ): List[PathEntry] = {
    //set orientation character
    val ori = if (alignment.isForward()) '+' else '-'
    //get alignment coordinates on reference genome
    val ref_coords = (alignment.rcoords._1, alignment.rcoords._2 + 1)
    //get local projected genes, in order
    val projected = cov_filter(ref_coords, local_genes.filterOverlaps(ref_coords).toList).map(_._2)
    //adjust orientation
    (if (alignment.isForward()) projected else projected.reverse).map(x => new PathEntry(x, ori))
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
                      ): ((Int, Int), List[((Int, Int), Int)]) => List[((Int, Int), Int)] = (read_range, genes) => {
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
