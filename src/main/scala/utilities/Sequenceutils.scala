package utilities

/**
  * Author: Alex N. Salazar
  * Created on 16-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object Sequenceutils {

  /**
    * Curried function to compute X-sequence content (i.e. GC or AT content). Requires a set of all desired
    * characters (i.e. Set('A','T') for AT-content). Case tolerable. Returns percentage in the given sequence.
    */
  def computeXcontent(xcontent: Set[Char]): String => Double = seq => {
    //iterate thorugh each nucleotide and compute x-sequence content
    val (content_count, total_count) = seq.foldLeft((0,0)){case((acc_content_count, acc_total_count),nt) => {
      ((if(xcontent(nt.toUpper)) (acc_content_count + 1) else acc_content_count), acc_total_count+1)
    }}
    //compute percentage
    content_count.toDouble / total_count
  }

}
