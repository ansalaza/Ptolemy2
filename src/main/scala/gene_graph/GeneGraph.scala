package gene_graph

import java.io.{File, PrintWriter}

import utilities.GeneGraphUtils.{GeneGraph, Paths, empty_gene_graph, empty_paths}
import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.GFFutils.parseMultiGFF2FingerTree
import utilities.PAFutils.{curateAlignmentsPerSeq}
import utilities.GFAutils.GFAwriter
import utilities.GeneProjectionUtils.{FingerTree}
import utilities.AlignmentGeneGraph.alignmentGeneGraph
import utilities.GenomeGeneGraph.ptolemyGeneGraph

/**
  * Author: Alex N. Salazar
  * Created on 7-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GeneGraph extends GFAwriter {

  case class Config(
                     genomeAlignments: File = null,
                     readAlignments: File = null,
                     gffFile: File = null,
                     outputDir: File = null,
                     minMapq: Int = 10,
                     minAlignmentCov: Double = 0.5,
                     featureTypes: String = null,
                     prefix: String = null,
                     nameTag: String = null,
                     excludeFile: File = null,
                     showWarning: Boolean = false)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("gene-graph") {
      opt[File]('g', "gff-files") required() action { (x, c) =>
        c.copy(gffFile = x)
      } text ("Tab-delimited file cantaining full path all GFF3-formatted files, one per line.")
      opt[String]("feature-types") required() action { (x, c) =>
        c.copy(featureTypes = x)
      } text ("Feature types in the GFF3 file to analyse as a string argument, comma-separated (i.e.the " +
        "string: 'gene,tRNA)")
      opt[String]("name-tag") required() action { (x, c) =>
        c.copy(nameTag = x)
      } text ("Name tag storing gene name description in the attributes column (i.e. the string: 'product=')")
      opt[File]('w', "wg-alignments") action { (x, c) =>
        c.copy(genomeAlignments = x)
      } text ("Whole genome alignments in PAF-format.")
      opt[File]('r', "read-alignments") action { (x, c) =>
        c.copy(readAlignments = x)
      } text ("Long-read alignments in PAF-format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      note("\nOPTIONAL\n")
      opt[Int]("min-mapq") action { (x, c) =>
        c.copy(minMapq = x)
      } text ("Process only alignments with at least this MAPQ (default is 10).")
      opt[Double]("alignment-coverage") action { (x, c) =>
        c.copy(minAlignmentCov = x)
      } text ("Mininum alignment coverage for a genes to be processed (default is 0.5).")

      opt[File]("exclude-sequence") action { (x, c) =>
        c.copy(excludeFile = x)
      } text ("Exclude sequence with the same IDs in the provided file (one per line).")
      opt[Unit]("show-warning") action { (x, c) =>
        c.copy(showWarning = true)
      } text ("Output warnings (turned-off by default).")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.gffFile)
      if (config.genomeAlignments == null && config.readAlignments == null) {
        assert(false, "Provide a genome alignmnet or read alignment file.")
      }
      if (config.genomeAlignments != null) verifyFile(config.genomeAlignments)
      if (config.readAlignments != null) verifyFile(config.readAlignments)
      constructGeneGraph(config)
    }
  }

  def constructGeneGraph(config: Config): Unit = {
    /* START: SET CURRIED METHODS/FUNCTIONS AND MISCELLANEOUS INFORMATION */
    //get feature types in a set
    val feature_types = config.featureTypes.replaceAll("\\s+", "").split(",").toSet
    println(timeStamp + "Found the following feature types: " + feature_types.mkString(","))
    //set curried multi-gff3 parser
    val parseGFFfiles = parseMultiGFF2FingerTree(feature_types, config.nameTag) _
    /* END: SET CURRIED METHODS/FUNCTIONS AND MISCELLANEOUS INFORMATION */

    /* START: PARSE GFF3-FORMATTED FILE INTO FINGER TREE DATA STRUCTURES */
    val (global_fingertree, global_description, global_origin) = {
      //first iterate through each gff and verify file is valid
      val gff_files = openFileWithIterator(config.gffFile).toList.map(file => {
        //create file
        val gff_file = new File(file)
        //verify
        verifyFile(gff_file)
        //return
        gff_file
      })
      parseGFFfiles(gff_files, config.showWarning)
    }
    println(timeStamp + "Found " + global_fingertree.size + " sequences with a total of " +
      global_description.size + " genes")
    /* END: PARSE GFF3-FORMATTED FILE INTO FINGER TREE DATA STRUCTURES */

    /**
      * Create gene graph from whole genome and long-read alignment. First create gene-graph from genome alignments
      * (note that an empty graph is returned if no whole genome alignments were provided). The resulting graph is
      * then updated with the long-read alignments (note not updates are made if no long-read alignments were provided).
      *
      * Returns 6-tuple: (final gene graph, genome paths (can be empty), node coverage from alignments (can be empty),
      * edge coverage from alignments (can be empty), temporary file contaning gene projecton individual reads (can
      * be null).
      */
    val (
      gene_graph,
      genome_paths,
      genome_node_coverage,
      genome_edge_coverage,
      alignment_node_coverage,
      alignment_edge_coverage,
      tmp_projection) = {
      //open and curate whole genome alignments, if any
      val (genome_gene_graph, genome_paths, sa_mappings, gnode_coverage, gedge_coverage) =
        createGenomeGeneGraph(global_fingertree, config)
      //open and curate long-read alignments, if any
      val g = createAlignmentGeneGraph(global_fingertree, config, genome_gene_graph, sa_mappings)
      //return final gene graph, genome paths,
      (g._1, genome_paths, gnode_coverage, gedge_coverage, g._2, g._3, g._4)
    }

    //total nodes
    val total_nodes = gene_graph.keys.size
    //total edges
    val total_edges = gene_graph.values.toList.map(_.size).sum
    println(timeStamp + "Constructed final gene graph of " + total_nodes + " nodes and " + total_edges + " edges")

    println(timeStamp + "Writing to disk")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/gene_graph.gfa")
    //output graph to gfa
    ptolemyGraph2GFA(pw, gene_graph, genome_paths, genome_node_coverage, genome_edge_coverage,
      alignment_node_coverage, alignment_edge_coverage, tmp_projection)
    pw.close()
    //create id file
    val pw2 = new PrintWriter(config.outputDir + "/gene_graph.node_ids.txt")
    global_description.toList.sortBy(_._1).foreach { case (node, (description, (start, end))) => {
      pw2.println(node + "\t" + global_origin(node) + "\t" + start + "\t" + end + "\t" + description)
    }
    }
    pw2.close
    println(timeStamp + "Successfully completed!")

  }

  /**
    * Method to create gene graphs from whole genome alignments. Will return empty graph and paths if no whole genome
    * alignments were provided
    *
    * @param fingertrees Finger tree data structures of all the genomes
    * @param config      configuration object
    * @return 5-tuple as (GeneGraph, Paths, SA-mappings, Node coverage, Edge coverage). SA-mappings is map(old id ->
    *         new id) for the new assigned id for each ortholog cluster
    */
  def createGenomeGeneGraph(fingertrees: Map[String, FingerTree],
                            config: Config): (GeneGraph, Paths, Map[Int, Int], Map[Int, Int], Map[(Int, Int), Int]) = {
    //user did not provide whole genome alignments
    if (config.genomeAlignments == null) {
      println(timeStamp + "Skipping gene graph construction from whole genome alignments")
      (empty_gene_graph, empty_paths, Map.empty[Int, Int], Map.empty[Int, Int], Map.empty[(Int,Int), Int])
    }
    //user provided whole genome alignments, construct gene graphs
    else {
      println(timeStamp + "Constructing native gene graphs from genomes")
      /* START: CURATE GENOME ALIGNMENTS*/
      println(timeStamp + "Curating whole genome alignments")
      //initialize id to genome map by parsing genome alignments file
      val (_id2genome, initial_id) = openFileWithIterator(config.genomeAlignments).foldLeft(List[(String, Int)]()) {
        case (genomes, line) => {
          //split line
          val split = line.split("\t")
          //add genome names and their sizes to list
          (split(5), split(6).toInt) :: ((split(0), split(1).toInt) :: genomes)
        }
      }
        //get all unique entries
        .toSet.toList
        //iterate and assigne unique IDs
        .foldLeft((List[(String, (Int, Int))](), 0)) { case ((genome_ids, id), (name, size)) => {
        ((name, (id, size)) :: genome_ids, id + 1)
      }
      }
      //curate whole genome alignments
      val (id2genome, curated_alignments) = {
        //curate alignments
        val tmp = curateAlignmentsPerSeq(config.genomeAlignments, config.minMapq, initial_id, _id2genome.toMap)
        //if no exclusive file found, move on
        if (config.excludeFile == null) tmp
        //exclude specified sequences
        else {
          //verify file
          verifyFile(config.excludeFile)
          //get all read ids to exclude
          val excluded_seqs = openFileWithIterator(config.excludeFile).toList.toSet
          println(timeStamp + "Found " + excluded_seqs.size + " sequences to exclude")
          //get corresponding ids for reads to exclude
          val excluded_ids = tmp._1.toList.foldLeft(List[Int]())((ids, read) => {
            //only add ids if they belong to the exclude set
            if (!excluded_seqs(read._2._1)) ids else read._1 :: ids
          }).toSet
          //filter out alignments
          (tmp._1.filterNot(x => excluded_seqs(x._2._1)), tmp._2.filterNot(x => excluded_ids(x._1)))
        }
      }
      /* END: CURATE GENOME ALIGNMENTS*/
      //construct ptolemy gene graph
      println(timeStamp + "Constructing gene graph from whole genome alignments")
      ptolemyGeneGraph(fingertrees, curated_alignments, config.minAlignmentCov, id2genome)
    }
  }

  def createAlignmentGeneGraph(fingertrees: Map[String, FingerTree],
                               config: Config,
                               genome_gene_graph: GeneGraph = empty_gene_graph,
                               sa_mappings: Map[Int, Int] = Map.empty[Int, Int]
                              ): (GeneGraph, Map[Int, Int], Map[(Int, Int), Int], File) = {
    //user did not provide read alignments
    if (config.readAlignments == null) {
      println(timeStamp + "Skipping gene graph construction from long-read alignments")
      (genome_gene_graph, Map(), Map(), null)
    }
    //user provided read alignments, create/update gene graph
    else {
      println(timeStamp + "Curating long-read alignments")
      //set temporary file to store projections of individual reads
      val tmp_file = new File(config.outputDir + "/.tmp_projections.gfa")
      //open and curate alignments
      val (id2readname, curated_alignments) = {
        //curate alignments
        val tmp = curateAlignmentsPerSeq(config.readAlignments, config.minMapq)
        //if no exclusive file found, move on
        if (config.excludeFile == null) tmp
        else {
          //verify file
          verifyFile(config.excludeFile)
          //get all read ids to exclude
          val excluded_reads = openFileWithIterator(config.excludeFile).toList.toSet
          println(timeStamp + "Found " + excluded_reads.size + " reads to exclude")
          //get corresponding ids for reads to exclude
          val excluded_ids = tmp._1.toList.foldLeft(List[Int]())((ids, read) => {
            //only add ids if they belong to the exclude set
            if (!excluded_reads(read._2._1)) ids else read._1 :: ids
          }).toSet
          //filter out alignments
          (tmp._1.filterNot(x => excluded_reads(x._2._1)), tmp._2.filterNot(x => excluded_ids(x._1)))
        }
      }
      println(timeStamp + "Curated alignments from a total of " + id2readname.size + " sequences")
      println(timeStamp + "Constructing gene graph")
      //construct gene graph
      val results = alignmentGeneGraph(fingertrees, curated_alignments, config.minAlignmentCov, id2readname,
        genome_gene_graph, sa_mappings, tmp_file)
      (results._1, results._2, results._3, tmp_file)
    }
  }
}
