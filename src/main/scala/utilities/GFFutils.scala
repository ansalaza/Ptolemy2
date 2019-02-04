package utilities

import java.io.File

import utilities.FileHandling.{openFileWithIterator, timeStamp}
import utilities.GeneProjectionUtils.{FingerTree, empty_fingertree}


/**
  * Author: Alex N. Salazar
  * Created on 13-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GFFutils {

  /**
    * Case class for storing a gene
    * @param id Unique gene ID
    * @param ori Orientation of gene
    */
  case class Gene(id: Int, ori: Char) {
    /**
      * Method to check whether give path entry is in forward orientation
      * @return Boolean
      */
    def isForward(): Boolean = ori == '+'

    /**
      * Reverse orientation of geme entry
      * @return PathEntry
      */
    def reverse(): Gene = new Gene(id, (if(isForward) '-' else '+'))
  }

  /**
    * Missing gene description
    */
  private val missing_gene_description = "null"

  /**
    * Case class for storing a GFF gene annotation as an object
    *
    * @param chrm
    * @param feature
    * @param start
    * @param end
    * @param ori
    * @param name
    */
  private case class GFFLine(chrm: String, feature: String, start: Int, end: Int, ori: Char, name: String)

  /**
    * Method to parse multiple GFF3-formatted files into a map of FingerTree data structures where the key is the
    * name of the contig/chromosome/scaffold. Curried function (first set of parameters are for initialization purposes)
    *
    * @param feature_types Set of strings representing which lines in the GFF3 file to process (e.g. 'CDS', 'gene')
    * @param name_tag      Tag in the description column indicating name/description of gene. (e.g. 'Name=')
    * @param gff_files     List of GFF3-formatted files as File objects
    * @param showWarning
    * @return Returns a 3-tuple:
    *         1 => Map(Sequence ID -> Finger tree)
    *         2 => Map(Gene ID -> (Gene description, (Start, End))
    *         3 => Map(Gene ID -> Sequence ID)
    */
  def parseMultiGFF2FingerTree(feature_types: Set[String], name_tag: String)
                              (gff_files: List[File],
                               showWarning: Boolean
                              ): (Map[String, FingerTree], Map[Int, (String, (Int, Int))], Map[Int, String]) = {
    //iterate through each gff and build global fingertree and gene identifier maps
    val (_global_fingertree, _global_description, _global_origin, _global_id) = {
      //iterate through each gff and construct local finger trees and identifiers
      gff_files.foldLeft((Map[String, FingerTree](), Map[Int, (String, (Int, Int))](), Map[Int, String](), 0)) {
        case ((fgt, description, genome_origin, id), gff) => {
          //get figertree, map, and sequence id for local gff
          val (updated_fgt, geneid2description, geneid2genome, last_id) =
            parseGFF2FingerTree(gff, feature_types, name_tag, showWarning, id, fgt)
          //add local gff's finger tree
          (updated_fgt,
            //add local gene identifiers
            geneid2description.toList.foldLeft(description.toList)((b, a) => a :: b).toMap,
            //add local gene ids to genome
            geneid2genome.toList.foldLeft(genome_origin.toList)((b, a) => a :: b).toMap,
            //Continue with recent added gene id
            last_id)
        }
      }
    }
    (_global_fingertree, _global_description, _global_origin)
  }

  /**
    * Method to parse GFF file and add annotations to a finger tree data structure. Assumes only one
    * sequence/contig/scaffold is present in the file.
    *
    * @param file        GFF3-formatted file
    * @param label       Set of labels correponding to which type of annotations to process
    * @param name_tag    Tag in the GFF3 attributes column storing gene name description
    * @param showWarning Output warnings about non-GFF3 formatted lines
    * @return A 5-tuple: (FingerTree data structure, Map(assigned id -> gene description), Map(assigned id -> genome),
    *         gene ID to cont., sequence ID)
    */

  def parseGFF2FingerTree(file: File,
                          label: Set[String],
                          name_tag: String,
                          showWarning: Boolean,
                          id: Int = 0,
                          inititial_fingertree: Map[String, FingerTree] = Map()
                         ): (Map[String, FingerTree], Map[Int, (String, (Int, Int))], Map[Int, String], Int) = {

    /**
      * Function to parse gff line and return GFFLine object
      *
      * @return GFFLine object
      */
    def toGFFLine: Array[String] => GFFLine = line => {
      //orientation of gene annotation
      val orientation = line(6).head
      assert(orientation == '+' || orientation == '-', "Unrecognized orientation: " + line)
      //start and end coordinates, 0-based
      val (start,end) = (line(3).toInt - 1, line(4).toInt)
      assert(start < end, "Unexpected start/end coordinates: " + line)
      //if there are less than 9 columns, non-GFF3 format
      new GFFLine(line(0), line(2), start, end, orientation, getGeneName(line(8)))
    }

    /**
      * Method to parse out the gene name in the attributes column of a gff3 file
      *
      * @param attribute Attribute column from a GFF3 formatted file
      * @return String
      */
    def getGeneName(attribute: String): String = {
      //split attribute line
      val annotation = attribute.split(";")
      //attempt to find name tag
      val gene_name = annotation.find(_.startsWith(name_tag))
      //return harcoded message if gene name description cannot be found
      if (gene_name.isEmpty) missing_gene_description
      //else return gene name description
      else gene_name.get.substring(gene_name.get.indexOf(name_tag) + name_tag.size)
    }

    //iterate through each line and add annotations to finger tree
    val (fingertree, id2interval, id2genome, gene_id) = openFileWithIterator(file)
      .foldLeft((inititial_fingertree, Map[Int, (String, (Int, Int))](), Map[Int, String](), id)) {
        //PASS ACC. MAP OF CHRM -> FINGERTREE
        case ((_fingertree, _geneid2description, _geneid2genome, _gene_id), line) => {
          //move on if its a comment line
          if (line.startsWith("#")) (_fingertree, _geneid2description, _geneid2genome, _gene_id)
          else {
            //count number of columns
            val columns = line.split("\t")
            //non-GFF3 formatted file
            if (columns.size < 9) {
              //log message
              if (showWarning) println(timeStamp + "----WARNING: found a non-GFF-formatted line. Skipping: " + line)
              //move on
              (_fingertree, _geneid2description, _geneid2genome, _gene_id)
            }
            //not a line corresponding to any of the labels specified
            else if (!label(columns(2))) (_fingertree, _geneid2description, _geneid2genome, _gene_id)
            //add to finger tree
            else {
              //parse gff line
              val gene = toGFFLine(columns)
              //fetch current finger tree of chrm and update
              val current = _fingertree.getOrElse(gene.chrm, empty_fingertree)
              //add gene to finger tree
              (_fingertree + (gene.chrm -> current.+((gene.start, gene.end), new Gene(_gene_id, gene.ori))),
                //add gene id to description
                _geneid2description + (_gene_id -> (gene.name, (gene.start, gene.end))),
                //add gene id to genome, increment gene id, return last sequence id added
                _geneid2genome + (_gene_id -> gene.chrm), _gene_id + 1)
            }
          }
        }
      }
    (fingertree, id2interval, id2genome, gene_id)
  }
}
