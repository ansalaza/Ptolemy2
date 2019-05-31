package utilities

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp}
import utilities.PAFutils.{PAFentry, toPAFentry}
import utilities.BlastUtils.{BLASTentry, toBlast}

/**
  * Author: Alex N. Salazar
  * Created on 20-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object MetaDataUtils {

  case class Description(id: Int, ref: String, start: Int, end: Int, ori: Char, description: String, geneID: String) {
    override def toString(): String = List(id, ref, start + 1, end, ori, description, geneID).mkString("\t")
  }

  /**
    * Function to load a file into  a list of Description objects
    *
    * @return List[Description]
    */
  def loadDescriptions: File => List[Description] = file => {
    /**
      * Function to compute line to Descriptions object
      *
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
    *
    * @return Map[Int, String]
    */
  def loadLabels: File => Map[Int, List[String]] = file => {
    //create graph as map of node ID -> list of labels
    openFileWithIterator(file).foldLeft(List[(Int, String)]())((acc, line) => {
      val columns = line.split("\t")
      (columns.head.toInt, columns(1)) :: acc
    }).groupBy(_._1).mapValues(_.map(_._2).distinct)
  }

  def loadAlignmentAsOrientations(file: File, isBlast: Boolean): Map[(String, String), Char] = {
    openFileWithIterator(file).toList.foldLeft(Map[(String, String), (Char, Int)]())((map, line) => {
      //parse alignment as either PAF or BLAST line
      val alignment: Either[PAFentry, BLASTentry] = if (isBlast) Right(toBlast(line)) else Left(toPAFentry(line))
      //get current alignment length
      val alignment_length = if(alignment.isLeft) alignment.left.get.align_length else alignment.right.get.alength
      //set alignment tuple of (query, subj)
      val atuple = {
        if (alignment.isLeft) (alignment.left.get.qname, alignment.left.get.rname)
        else (alignment.right.get.qname, alignment.right.get.rname)
      }
      //move on if self-alignment
      if (atuple._1 == atuple._2) map
      else {
        //set alignment orientation
        val ori = if (alignment.isLeft) alignment.left.get.ori else alignment.right.get.ori
          //assert(!(map.contains(atuple) || map.contains(atuple)), "The two sequences are already in the map: " +
            //atuple + ". Based on alignment: " + line)
        //add to map
        if(!map.contains(atuple)) map + (atuple -> (ori, alignment_length))
        //existing multiple alignments to same gene
        else {
          //get new ori and legnth
          val (new_ori, new_length) = map(atuple)
          if(new_ori != ori) println(timeStamp + "WARNING: dual alignments with different orientations: " + line)
          //move on
          if(alignment_length >= new_length) map
          //update with updated orientation
          else map + (atuple -> (new_ori, new_length))
        }
      }
    }).mapValues(_._1)
  }

  def loadGeneFamilies: File => Map[Int, List[Int]] = file => {
    openFileWithIterator(file).foldLeft(List[(Int, Int)]())((acc, line) => {
      val columns = line.split("\t")
      assert(columns.size == 2, "Found multiple columns in line: " + columns.toList)
      (columns.head.toInt, columns(1).toInt) :: acc
    }).groupBy(_._2).mapValues(_.map(_._1).toList.distinct)
  }

}
