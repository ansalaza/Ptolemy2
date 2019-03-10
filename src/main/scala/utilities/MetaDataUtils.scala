package utilities

import java.io.File
import utilities.FileHandling.openFileWithIterator

/**
  * Author: Alex N. Salazar
  * Created on 20-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object MetaDataUtils {

  case class Description(id: Int, ref: String, start: Int, end: Int, ori: Char, description: String, geneID: String) {
    override def toString(): String = List(id, ref, start, end, ori, description, geneID).mkString("\t")
  }

  /**
    * Function to load a file into  a list of Description objects
    * @return List[Description]
    */
  def loadDescriptions: File => List[Description] = file => {
    /**
      * Function to compute line to Descriptions object
      * @return Description
      */
    def toDescription: String => Description = line => {
      //get columns
      val columns = line.split("\t")
      new Description(columns(0).toInt, columns(1), columns(2).toInt, columns(3).toInt, columns(4).head, columns(5),
        columns(6))
    }
    //iterate through lines and create Description instances
    openFileWithIterator(file).toList.map(toDescription(_))
  }

  /**
    * Function to load a labels file containing node id to label
    * @return Map[Int, String]
    */
  def loadLabels: File => Map[Int, List[String]] = file => {
    //create graph as map of node ID -> list of labels
    openFileWithIterator(file).foldLeft(List[(Int,String)]())((acc, line) => {
      val columns = line.split("\t")
      (columns.head.toInt, columns(1)) :: acc
    }).groupBy(_._1).mapValues(_.map(_._2).distinct)
  }

  def loadGeneFamilies: File => Map[Int, List[Int]] = file => {
    openFileWithIterator(file).foldLeft(List[(Int,Int)]())((acc, line) => {
      val columns = line.split("\t")
      assert(columns.size == 2, "Found multiple columns in line: " + columns.toList)
      (columns.head.toInt, columns(1).toInt) :: acc
    }).groupBy(_._2).mapValues(_.map(_._1).toList.distinct)
  }

}
