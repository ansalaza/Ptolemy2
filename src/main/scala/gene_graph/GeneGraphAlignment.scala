package gene_graph

import java.io.{File, PrintWriter}

import utilities.GeneGraphUtils.{GeneGraph, empty_gene_graph}
import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.GFFutils.parseMultiGFF2FingerTree
import utilities.PAFutils.curateAlignmentsPerSeq
import utilities.GFAutils.GeneGraphWriter
import utilities.GeneProjectionUtils.FingerTree
import utilities.AlignmentGeneGraph.alignmentGeneGraph

/**
  * Author: Alex N. Salazar
  * Created on 7-12-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GeneGraphAlignment extends GeneGraphWriter {

  case class Config(
                     readAlignments: File = null,
                     gffFile: File = null,
                     outputDir: File = null,
                     minMapq: Int = 10,
                     minMultiMap: Int = 100,
                     minDist: Int = 100,
                     minAlignmentCov: Double = 0.5,
                     featureType: String = null,
                     prefix: String = null,
                     nameTag: String = null,
                     idTag: String = "protein_id=",
                     excludeSeqFile: File = null,
                     excludeRegionFile: File = null,
                     splitOri: Boolean = false,
                     showWarning: Boolean = false)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("gene-graph") {
      opt[File]('g', "gff-files") required() action { (x, c) =>
        c.copy(gffFile = x)
      } text ("File containing full path to all corresponding GFF3-formatted files, one per line.")
      opt[File]('r', "read-alignments") required() action { (x, c) =>
        c.copy(readAlignments = x)
      } text ("Long-read alignments in PAF-format.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory. If it doesn't not exist, the directory will be created.")
      opt[String]("prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output file.")
      note("\nOPTIONAL")
      note("\nINPUT")
      opt[String]("feature-type") action { (x, c) =>
        c.copy(featureType = x)
      } text ("Feature types in the GFF3 file to analyse (default is 'CDS').")
      opt[String]("name-tag") action { (x, c) =>
        c.copy(nameTag = x)
      } text ("Name tag storing gene name description in the attributes column (default is 'product=')")
      opt[String]("id-tag") action { (x, c) =>
        c.copy(idTag = x)
      } text ("Name tag storing gene ID (default is 'protein_id=')")
      opt[File]("exclude-sequence") action { (x, c) =>
        c.copy(excludeSeqFile = x)
      } text ("Exclude sequence with the same IDs in the provided alignment(s) file (one per line).")
      opt[File]("exclude-regions") action { (x, c) =>
        c.copy(excludeRegionFile = x)
      } text ("Exclude reference (sub-)regions (one per line).")
      note("\nALGORITHM")
      opt[Unit]("split-ori") action { (x,c) =>
        c.copy(splitOri = true)
      } text ("Do not collapse homologous genes unless they are in the same orientation.")
      opt[Int]("min-mapq") action { (x, c) =>
        c.copy(minMapq = x)
      } text ("Process only alignments with at least this MAPQ (default is 10).")
      opt[Int]("min-overlap") action { (x, c) =>
        c.copy(minMultiMap = x)
      } text ("Minimum overlap length for an alignment to be considered multi-mapping (default is 100).")
      opt[Int]("min-dist") action { (x, c) =>
        c.copy(minDist = x)
      } text ("Minimum distance for 1D clustering of multimapping breakpoints (default is 100).")
      opt[Double]("alignment-coverage") action { (x, c) =>
        c.copy(minAlignmentCov = x)
      } text ("Mininum alignment coverage for a genes to be processed (default is 0.75).")
      opt[Unit]("show-warning") action { (x, c) =>
        c.copy(showWarning = true)
      } text ("Output warnings (turned-off by default).")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.gffFile)
      if (config.readAlignments != null) verifyFile(config.readAlignments)
      constructGeneGraph(config)
    }
  }

  def constructGeneGraph(config: Config): Unit = {
    /* START: SET CURRIED METHODS/FUNCTIONS AND MISCELLANEOUS INFORMATION */
    //set curried multi-gff3 parser
    val parseGFFfiles = parseMultiGFF2FingerTree(Set(config.featureType), config.nameTag, config.idTag) _
    /* END: SET CURRIED METHODS/FUNCTIONS AND MISCELLANEOUS INFORMATION */

    /* START: PARSE GFF3-FORMATTED FILE INTO FINGER TREE DATA STRUCTURES */
    val (global_fingertree, global_description) = {
      //first iterate through each gff and verify file is valid
      val gff_files = openFileWithIterator(config.gffFile).toList.map(line => {
        //make file
        val gff_file = new File(line)
        //verify
        verifyFile(gff_file)
        //return
        gff_file
      })
      parseGFFfiles(gff_files, config.showWarning)
    }
    println(timeStamp + "Found " + global_fingertree.size + " sequences with a total of " + global_description.size +
      " genes")
    /* END: PARSE GFF3-FORMATTED FILE INTO FINGER TREE DATA STRUCTURES */

    /**
      * Create gene graph from long-read alignment.
      * Returns 6-tuple: (final gene graph, genome paths (can be empty), node coverage from alignments (can be empty),
      * edge coverage from alignments (can be empty), temporary file contaning gene projecton individual reads (can
      * be null).
      */
    val (gene_graph, alignment_node_coverage, alignment_edge_coverage, tmp_projection) = {
      //open and curate long-read alignments, if any
      val g = createAlignmentGeneGraph(global_fingertree, config)
      //return final gene graph, genome paths,
      (g._1, g._2, g._3, g._4)
    }

    //total nodes
    val total_nodes = gene_graph.keys.size
    //total edges
    val total_edges = gene_graph.values.toList.map(_.size).sum
    println(timeStamp + "Constructed final gene graph of " + total_nodes + " nodes and " + total_edges + " edges")

    println(timeStamp + "Writing to disk")
    //create output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".gfa")
    //output graph to gfa
    geneGraph2GFA(pw, gene_graph, None, alignment_node_coverage, alignment_edge_coverage, tmp_projection)
    pw.close()
    //create id file
    val pw2 = new PrintWriter(config.outputDir + "/" + config.prefix + ".node_ids.txt")
    global_description.toList.sortBy(_._1).foreach(x => pw2.println(x._2.toString()))
    pw2.close
    println(timeStamp + "Successfully completed!")

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
        val tmp = curateAlignmentsPerSeq(config.readAlignments, config.minMapq, config.minMultiMap, config.minDist)
        //if no exclusive file found, move on
        if (config.excludeSeqFile == null) tmp
        else {
          //verify file
          verifyFile(config.excludeSeqFile)
          //get all read ids to exclude
          val excluded_reads = openFileWithIterator(config.excludeSeqFile).toList.toSet
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
      /**
      //get all multimap alignment intervals
      val low_qual_regions = curated_alignments.foldLeft(Map[String, List[(Int,Int)]]()){
        //iterate through curated alignments
        case (multimaps, (id, alignment_intervals)) => {
          //for each curated alignment interval, check number of alignments
          alignment_intervals.foldLeft(multimaps)((acc, alignments) =>
            //only add alignments if they are a multimapped region
            if(alignments.size == 1) acc
            else {
              //update map with current ref multimap coordinates
              alignments.foldLeft(acc)((a, alignment) => {
                println(alignment)
                val fetch = a.getOrElse(alignment.ref, List[(Int,Int)]())
                a + (alignment.ref -> (alignment.rcoords :: fetch))
              })
            }
          )
        }}//obtain canonical multimap regions
        .mapValues(x => longestOverlappingIntervals(x).sortBy(identity))
      //log low quality regions
      if(low_qual_regions.nonEmpty){
        println(timeStamp + "Found " + low_qual_regions.values.map(_.size).sum + " ambiguous regions. Omitting:")
        low_qual_regions.foreach(x => {
          println(timeStamp + "--" + x._1)
          x._2.foreach(y => println(timeStamp + "----" + y))
        })
      }
        */
      /**
      //update SA mapping , if needed
      val updated_sa_mapping = {
        //get all multimap alignment intervals
        val multimapped = curated_alignments.foldLeft(List[(String, List[List[Alignment]])]()){
          //iterate through curated alignments
          case (multimaps, (id, alignment_intervals)) => {
            //for each curated alignment interval, check number of alignments
          alignment_intervals.foldLeft(multimaps)((acc, alignments) =>
            //only add alignments if they are a multimapped region
            if(alignments.size == 1) acc
            else pairwiseTargetAlignments(alignments).foldLeft(acc)((a, b) => b :: a)
            )
        }}
        if(multimapped.isEmpty) sa_mappings
        else {
          println(timeStamp + "Found " + multimapped.size + " multimapping regions")
          println(timeStamp + "Updating SA mapping")
          //create multimapped-based sa graph
          val multimap_based_sa = {
            //create ref name -> local ID
            val ref2ID = multimapped.map(_._1).zipWithIndex.toMap
            //find syntenica anchors
            val sa = findSytenicAnchors(fingertrees, ref2ID.map(x => (x._2, (x._1, -1))), config.minAlignmentCov, false,
              multimapped.map(x => (ref2ID(x._1), x._2)))
            //create SA graph
            createSAgraph(sa)
          }
          //get local sa mappings
          val local_mappings = defineCanonicalIDs(getConnectedComponents(multimap_based_sa))
          println("LOCAL")
          println(local_mappings)
          println("PREVIOUS")
          println(sa_mappings)
          println(timeStamp + "Defined " + local_mappings.size + " new canonical mappings")
          local_mappings.foldLeft(sa_mappings)((map, local) => map + (local._1 -> local._2))

        }}
      println("UPDATED")
      updated_sa_mapping.foreach(println)
        */
      println(timeStamp + "Constructing gene graph")
      //construct gene graph
      val results = alignmentGeneGraph(fingertrees, curated_alignments, config.minAlignmentCov, id2readname,
        genome_gene_graph, sa_mappings, tmp_file)
      (results._1, results._2, results._3, tmp_file)
    }
  }
}




/*

  Archived code:

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
        val tmp = curateAlignmentsPerSeq(config.genomeAlignments, config.minMapq, config.minMultiMap, config.minDist,
          initial_id, _id2genome.toMap)
        //if no exclusive file found, move on
        if (config.excludeSeqFile == null) tmp
        //exclude specified sequences
        else {
          //verify file
          verifyFile(config.excludeSeqFile)
          //get all read ids to exclude
          val excluded_seqs = openFileWithIterator(config.excludeSeqFile).toList.toSet
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
      //identify syntenic anchors
      val sa_graph = {
        //create SA graph
        createSAgraph(
          //identify all syntenic anchors
        findSytenicAnchors(fingertrees, id2genome, config.minAlignmentCov, config.splitOri,
          curated_alignments)
        )
      }
      //construct ptolemy gene graph
      println(timeStamp + "Constructing gene graph from whole genome alignments")
      ptolemyGeneGraph(fingertrees, sa_graph)
    }
  }

  */