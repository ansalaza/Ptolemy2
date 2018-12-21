package cli

import de.sciss.fingertree.RangedSeq

/**
  * Author: Alex N. Salazar
  * Created on 21-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object play {

  def main(args: Array[String]): Unit = {
    val sq = RangedSeq(
      (1685, 1750) -> "Bach",
      (1866, 1925) -> "Satie",
      (1883, 1947) -> "Russolo",
      (1883, 1965) -> "Varèse",
      (1910, 1995) -> "Schaeffer",
      (1912, 1992) -> "Cage"
    )(_._1, Ordering.Int)

    implicit class Names(it: Iterator[(_, _)]) {
      def names = it.map(_._2).mkString(", ")
    }

    println(sq.filterOverlaps((1700,1900)).toList)

    println(List(List(1,2,3), List(4,5), List(6)).flatten)

  }

}
