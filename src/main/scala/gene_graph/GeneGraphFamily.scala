package gene_graph

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile}
import utilities.GFFutils.Gene
import utilities.GeneGraphUtils.{createSAgraph, defineCanonicalIDs, getConnectedComponents, _}
import utilities.PAFutils.{PAFentry, toPAFentry}
import utilities.MetaDataUtils.{loadDescriptions, loadAlignmentAsOrientations}
import utilities.BlastUtils.{BLASTentry, toBlast}
import utilities.GFAutils.GeneGraphWriter

/**
  * Author: Alex N. Salazar
  * Created on 16-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GeneGraphFamily extends GeneGraphWriter {

  case class Config(
                     alignments: File = null,
                     descriptionsFile: File = null,
                     outputDir: File = null,
                     dnaAlignments: File = null,
                     minCov: Double = 0.75,
                     minEval: Double = 0.000001,
                     minBitscore: Double = 50.0,
                     splitOri: Boolean = false,
                     prefix: String = null,
                     localID: Boolean = false,
                     isProtein: Boolean = false,
                     isBlast: Boolean = false,
                     dnaIsBlast: Boolean = false,
                     dnaIsLocal: Boolean = false,
                     metrics: Boolean = true)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("extract-genes") {
      opt[File]('a', "alignments") required() action { (x, c) =>
        c.copy(alignments = x)
      } text ("Pairwise alignments of DNA/protein gene sequences. Default format is DNA alignments in PAF-format. See" +
        " below parameters for protein-based and BLAST-formatted alignments.")
      opt[File]('d', "descriptions") required() action { (x, c) =>
        c.copy(descriptionsFile = x)
      } text ("Ptolemy's descriptions file (see README for more details).")
      opt[String]('p', "prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output files.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      note("\nOPTIONAL")
      note("\nINPUT")
      opt[Unit]("blast") action { (x, c) =>
        c.copy(isBlast = true)
      } text ("Alignments are formatted in BLAST-tabular format. Assumes BLAST was ran with the following " +
        "output parameter: '-outfmt \"6 qseqid qlen sseqid slen sstart send length pident evalue bitscore\"'")
      opt[Unit]("protein") action { (x, c) =>
        c.copy(isProtein = true)
      } text ("Protein-based alignments. Note that protein-based alignments are only compatible with BLAST-tabular " +
        "format (see '--blast' parameter).'")
      opt[Unit]("local-id") action { (x, c) =>
        c.copy(localID = true)
      } text ("Subject and target names input alignments are local gene/protein IDs (must match those in the " +
        "Ptolemy's descriptions file).")
      opt[File]("dna-alignments") action { (x, c) =>
        c.copy(dnaAlignments = x, splitOri = true)
      } text ("Only compatible with protein-based alignments. Provide DNA-based pairwise gene alignments  to evaluate" +
        " alignments via the '--split-ori' parameter (see below). PAF-format by default, change to BLAST via " +
        "'--dna-blast' parameter (see below).")
      opt[Unit]("dna-blast") action { (x, c) =>
        c.copy(dnaIsBlast = true)
      } text ("DNA-alignments are in same format specified in '--blast' parameter.")
      opt[Unit]("dna-local") action { (x, c) =>
        c.copy(dnaIsLocal = true)
      } text ("DNA-alignments use local protein/gene IDs (see '--local-id' parameter below).")
      note("\nHEURISTICS")
      opt[Unit]("split-ori") action { (x, c) =>
        c.copy(splitOri = true)
      } text ("Do not collapse homologous genes unless they are in the same orientation. Automatically determined " +
        "when using DNA-based alignments. Requires '--orientation' parameter when using protein-based alignments.")
      opt[Double]("min-cov") action { (x, c) =>
        c.copy(minCov = x)
      } text ("Minimum alignment coverage (default is 0.75).")
      opt[Double]("min-eval") action { (x, c) =>
        c.copy(minEval = x)
      } text ("Minimum alignment e-value (default is 0.000001). Only applicable when using BLAST-format.")
      opt[Double]("min-bitscore") action { (x, c) =>
        c.copy(minBitscore = x)
      } text ("Minimum alignment bitscore (default is 50.0). Only applicable when using BLAST-format")
    }
    parser.parse(args, Config()).map { config =>
      //check whether alignment file and output directory exist.
      verifyDirectory(config.outputDir)
      verifyFile(config.alignments)
      verifyFile(config.descriptionsFile)
      //sanity check for protein-based alignments
      if (config.isProtein) {
        if (config.splitOri) assert(config.dnaAlignments != null, "Using " +
          "protein-based alignments with '--split-ori' parameter turned on. Requires '--dna-alignments'.")
      }
      if(config.dnaAlignments != null) verifyFile(config.dnaAlignments)
      geneGraphFamily(config)
    }
  }

  def geneGraphFamily(config: Config): Unit = {
    //load descriptions file into id -> description
    val descriptions = loadDescriptions(config.descriptionsFile).map(x => x.id -> x).toMap
    //set local ID -> ptolemy id, if specified
    val l2p = if (!config.localID) Map[String, Int]() else descriptions.map(x => (x._2.geneID, x._1))
    //set orientations mapping, if provided
    val orientations = {
      if (!(config.splitOri && config.isProtein)) Map[(String, String), Char]()
      else {
        val tmp = loadAlignmentAsOrientations(config.dnaAlignments, config.dnaIsBlast)
        if(!config.dnaIsLocal) tmp else tmp.map(x => (l2p(x._1._1).toString, l2p(x._1._2).toString) -> x._2)
      }
    }.filter(x => (descriptions(x._1.toString).ref == descriptions(x._2.toString).ref))

    /**
      * Function to parse an alignment line (either PAF or BLAST-tabular)
      *
      * @return Option[(query name, ref name, orientation)]
      */
    def parseAlignment: String => Option[(String, String, Char)] = line => {
      /**
        * Function to parse a given line representing a PAF or BLAST-formatted line and returns either object with left
        * and right, respectively
        *
        * @return Either[PAFentry,BLASTentry]
        */
      def parseLine(): Either[PAFentry, BLASTentry] = if (config.isBlast) Right(toBlast(line)) else Left(toPAFentry(line))

      /**
        * FUnction to compute the alignment coverage of either PAFentry or BLASTentry
        *
        * @return Double
        */
      def computeCov: Either[PAFentry, BLASTentry] => Double = a => {
        val (x, y) = {
          if (a.isLeft) (a.left.get.query_length, a.left.get.ref_length)
          else (a.right.get.qlength, a.right.get.rlength)
        }
        min(x / y.toDouble, y / x.toDouble)
      }

      /**
        * Load either a PAFentry or BLASTentry as an option 3-tuple: (query name, ref name, orientation)
        *
        * @return Option[(String,String,Char)]
        */
      def loadEither: Either[PAFentry, BLASTentry] => Option[(String, String, Char)] = either => {
        if (either.isLeft) Option(either.left.get.qname, either.left.get.rname, either.left.get.ori)
        else {
          //set alignment orientation for current alignment
          val true_orientation = {
            //DNA-based alignment or protein-based alignment but no split ori
            if (!config.isProtein || !config.splitOri) either.right.get.ori
            //protein-based alignmetn with split ori
            else {
              //set alignment tuple
              val tuple = {
                val ids = (either.right.get.qname, either.right.get.rname)
                println(ids)
                if(!config.localID) ids else (l2p(ids._1).toString, l2p(ids._2).toString)
              }
              //get alignment orientation
              orientations.getOrElse(tuple, orientations(tuple.swap))
            }
          }
          //return option 3-tuple
          Option(either.right.get.qname, either.right.get.rname, true_orientation)
        }
      }

      /**
        * Function to determine whether a given alignment either PAF or BLAST-formatted is valid (passes heuristics)
        *
        * @return Boolean
        */
      def isValid: Either[PAFentry, BLASTentry] => Boolean = either => {
        //compute coverage
        val cov = computeCov(either)
        //PAF or BLAST of DNA-based alignment
        if (either.isLeft || !config.isProtein) cov > config.minCov
        //protein based alignment
        else {
          //load as right
          val right = either.right.get
          //check heuristics
          (right.eval < config.minEval && right.bitscore > config.minBitscore && cov > config.minCov)
        }
      }

      //parse alignment line
      val alignment = parseLine()
      //check whether is passes heuristic thresholds, if so return optional 3-tuple
      if (!isValid(alignment)) None else loadEither(alignment)
    }

    /**
      * Function to create a Gene object given a local/ptolemy ID
      *
      * @return Gene
      */
    def id2Gene: String => Gene = _id => {
      val id = if (!config.localID) _id.toInt else l2p(_id)
      new Gene(id, descriptions(id).ori)
    }

    /**
      * Function to determine whether a two given genes are in the same orientation based on their
      * DNA-based alignment orientation
      * @return Boolean
      */
    def isSameOri: (Gene, Gene, Char) => Boolean = (g1, g2, alignment_ori) => {
      (alignment_ori == '+' && g1.ori == g2.ori) || (alignment_ori == '-' && g1.ori != g2.ori)
    }

    println(timeStamp + "Loaded " + descriptions.size + " gene descriptions")
    println(timeStamp + "Clustering sequences based on alignments")
    //iterate through alignments and construct SA graph
    val sa = openFileWithIterator(config.alignments).foldLeft(List[(Gene, Gene)]())((sa, line) => {
      //parse alignment line
      val align = parseAlignment(line)
      //add to list only if alignment coverage suffice
      if (align.isEmpty) sa
      else {
        //set genes
        val (g1, g2) = (id2Gene(align.get._1), id2Gene(align.get._2))
        //genes can't be from the same sequence
        if (descriptions(g1.id).ref == descriptions(g2.id).ref) sa
        //split ori is not turned on, add to SA
        else if(!config.splitOri) (g1,g2) :: sa
        //split ori is turned on, determine orientation compatability
        else {
          //only add if genes are on the same orientation
          if(!isSameOri(g1, g2, align.get._3)) sa else (g1,g2) :: sa
        }
      }
    })
    println(timeStamp + "Building homology graph")
    //build graph
    val sa_graph = createSAgraph(sa)
    println(timeStamp + "Finding connected components")
    //get connected components
    val ccs = getConnectedComponents(sa_graph)
    //output metrics, if speciified
    if (config.metrics) {
      val pw_m = new PrintWriter(config.outputDir + "/" + config.prefix + ".homology_graph_connectivity.txt")
      pw_m.println("Ndegree\tCCsize")
      //compute ccs connectevitiy metrics
      ccs.foreach(cc => {
        //set cc size
        val size = cc.size.toDouble
        cc.foreach(node => {
          val edges = sa_graph.getOrElse(node, List[Int]())
          edges.forall(x => sa_graph(x).exists(_ == x))
          pw_m.println(edges.size / size + "\t" + size)
        })
      })
      pw_m.close
    }
    println(timeStamp + "--Found " + ccs.size + " connected components with an average size of " +
      (ccs.map(_.size).sum.toDouble / ccs.size))
    println(timeStamp + "Defining gene families")
    //set gene families
    val gene_families = defineCanonicalIDs(ccs)
    //set output file
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".gene_families.txt")
    //iterate through each gene family mapping and output
    gene_families.foreach { case (x, y) => pw.println(x + "\t" + y) }
    pw.close
    println(timeStamp + "Constructing genome-architectures")
    //reload gene description and convert to list of paths
    val genome_paths = {
      descriptions.values.groupBy(_.ref).mapValues(_.toList.sortBy(x => (x.start, x.end))
        .map(x => new Gene(gene_families.getOrElse(x.id, x.id), x.ori)))
    }.mapValues(x => (x, -1))
    println(timeStamp + "--Constructed architectures for " + genome_paths.size + " sequences")
    println(timeStamp + "Constructing gene-graph with coverage statistics")
    //construct gene graph along with node and edge coverage
    val (graph, ncov, ecov) = {
      val tmp = genome_paths.foldLeft((empty_gene_graph, empty_node_coverage, empty_edge_coverage)) {
        case ((_graph, _ncov, _ecov), (seq, (path, size))) => updateGeneGraph(path, false, _graph, _ncov, _ecov)
      }
      (tmp._1.mapValues(_.distinct), tmp._2, tmp._3)
    }
    println(timeStamp + "--Constructed gene graph with " + computeTotalNodes(graph) + " nodes and " + (ecov.size / 2) +
      " distinct edges")
    println(timeStamp + "Writing to disk")
    val pw2 = new PrintWriter(config.outputDir + "/" + config.prefix + ".gfa")
    geneGraph2GFA(pw2, graph, Option(genome_paths), ncov, ecov)
    pw2.close
    println(timeStamp + "Successfully completed!")

  }

  /**
    * Local function to compute min for double's
    *
    * @return
    */
  def min: (Double, Double) => Double = (x, y) => if (x < y) x else y


}
