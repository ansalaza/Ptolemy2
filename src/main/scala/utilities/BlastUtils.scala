package utilities

/**
  * Author: Alex N. Salazar
  * Created on 20-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object BlastUtils {

  /**
    * Case class for storing a blast alignment. Assumes alignment came from following format:
    * -outfmt "6 qseqid qlen sseqid slen length pident evalue bitscore"
    *
    * @param qname
    * @param qlength
    * @param rname
    * @param rlength
    * @param identity
    * @param eval
    * @param bitscore
    */
  case class BLASTentry(qname: String, qlength: Int, rname: String, rlength: Int, alength: Int, identity: Double,
  eval: Double, bitscore: Double)

  def toBlast: String => BLASTentry = line => {
    val columns = line.split("\t")
    new BLASTentry(columns.head, columns(1).toInt, columns(2), columns(3).toInt, columns(4).toInt,
      columns(5).toDouble, columns(6).toDouble, columns(7).toDouble)
  }

}
