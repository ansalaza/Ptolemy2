package gene_graph

import java.io.{File, PrintWriter}

import utilities.FileHandling.{openFileWithIterator, timeStamp, verifyDirectory, verifyFile, getFileExtension}
import utilities.GFFutils.Gene
import utilities.GeneGraphUtils.{createSAgraph, defineCanonicalIDs, getConnectedComponents, _}
import utilities.PAFutils.toPAFentry
import utilities.MetaDataUtils.loadDescriptions
import utilities.BlastUtils.toBlast
/**
  * Author: Alex N. Salazar
  * Created on 16-2-2019
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object GeneGraphFamily {

  case class Config(
                     alignments: File = null,
                     descriptionsFile: File = null,
                     outputDir: File = null,
                     minCov: Double = 0.75,
                     minEval: Double = 0.000001,
                     minBitscore: Double = 50.0,
                     splitOri: Boolean = false,
                     prefix: String = null,
                     localID: Boolean = false,
                     isBlast: Boolean = false,
                     metrics: Boolean = false)

  def main(args: Array[String]) {
    val parser = new scopt.OptionParser[Config]("extract-genes") {
      opt[File]('a', "alignments") required() action { (x, c) =>
        c.copy(alignments = x)
      } text ("Pairwise alignments of gene sequences in PAF-format. For BLAST-alignments, see '--blast' parameter.")
      opt[File]('d', "descriptions") required() action { (x, c) =>
        c.copy(descriptionsFile = x)
      } text ("Descriptions file (see README).")
      opt[String]('p', "prefix") required() action { (x, c) =>
        c.copy(prefix = x)
      } text ("Prefix for output files.")
      opt[File]('o', "output-directory") required() action { (x, c) =>
        c.copy(outputDir = x)
      } text ("Output directory to store indexed genomes.")
      note("\nOPTIONAL")
      note("\nDNA-BASED")
      opt[Double]("min-cov") action { (x, c) =>
        c.copy(minCov = x)
      } text ("Minimum alignment coverage (default is 0.75).")
      opt[Unit]("split-ori") action { (x, c) =>
        c.copy(splitOri = true)
      } text ("Do not collapse homologous genes unless they are in the same orientation.")
      note("\nPROTEIN-BASED")
      opt[Unit]("protein-blast") action { (x,c) =>
        c.copy(isBlast = true)
      } text ("(Protein) alignments are in BLAST-tabular format. Filters alignments based on '--min-eval' and " +
        "'--bitscore'. Assumes output was formatted with the following command: " +
        "'-outfmt \"6 qseqid qlen sseqid slen length pident evalue bitscore\"'")
      opt[Double]("min-eval") action {(x,c) =>
        c.copy(minEval = x)
      } text ("Minimum alignment e-value (default is 0.000001)")
      opt[Double]("min-bitscore") action {(x,c) =>
        c.copy(minBitscore = x)
      } text("Minimum alignment bitscore (default is 50.0)")
      opt[Double]("min-cov") action { (x, c) =>
        c.copy(minCov = x)
      } text ("Minimum alignment coverage (default is 0.75).")
      note("\nINPUT")
      opt[Unit]("local-id") action { (x,c) =>
        c.copy(localID = true)
      } text ("Subject and target names in alignment file are local gene/protein IDs (must match those in the " +
        "Descriptions file).")
      note("\nOUTPUT")
      opt[Unit]("metrics") action { (x, c) =>
        c.copy(metrics = true)
      } text ("Write temporary metrics to disk.")
    }
    parser.parse(args, Config()).map { config =>
      //check whether output directory exists. If not, create it.
      verifyDirectory(config.outputDir)
      verifyFile(config.alignments)
      geneGraphFamily(config)
    }
  }

  def geneGraphFamily(config: Config): Unit = {
    /**
      * Function to parse an alignment line (either PAF or BLAST-tabular
      * @return  Option[(query name, ref name, orientation)]
      */
    def parseAlignment: String => Option[(String, String, Char)] = line => {
      def computeCov: (Int,Int) => Double = (x,y) => min(x / y.toDouble, y / x.toDouble)
      if(!config.isBlast){
        //parse PAF line
        val paf = toPAFentry(line)
        //get min alignment coverage
        val cov = computeCov(paf.query_length, paf.ref_length)
        //check coverage threshold
        if(cov < config.minCov) None else Option(paf.qname, paf.rname, paf.ori)
      } else {
        //parse blast line
        val blast = toBlast(line)
        //get min alignment coverage
        val cov = computeCov(blast.qlength, blast.rlength)
        //check eval and bitscore thresholds
        if(blast.eval > config.minEval || blast.bitscore < config.minBitscore || cov < config.minCov) None
        else Some(blast.qname, blast.rname, '+')
      }
    }
    //load descriptions file into id -> description
    val descriptions = loadDescriptions(config.descriptionsFile).map(x => x.id -> x).toMap
    //set local ID -> id, if specified
    val p2i = if(!config.localID) Map[String, Int]() else descriptions.map(x => (x._2.geneID, x._1))
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
        val (g1, g2) = {
          //set gene IDs
          val (id1, id2) =
            if(!config.localID) (align.get._1.toInt, align.get._2.toInt) else (p2i(align.get._1), p2i(align.get._2))
          //set genes
          (new Gene(id1, descriptions(id1).ori), new Gene(id2, descriptions(id2).ori))
        }
        //genes can't be from the same sequence
        if (descriptions(g1.id).ref == descriptions(g2.id).ref) sa
        else {
          //if orientation parameter is turned on and not a blast alignment
          if (!config.isBlast && config.splitOri) {
            //check orientations
            if ((align.get._3 == '+' && g1.ori == g2.ori) || (align.get._3 == '-' && g1.ori != g2.ori)) (g1, g2) :: sa
            else sa
          }
          //add to list
          else (g1, g2) :: sa
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
    if(config.metrics){
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
    val pw = new PrintWriter(config.outputDir + "/" + config.prefix + ".txt")
    //iterate through each gene family mapping and output
    gene_families.foreach { case (x, y) => pw.println(x + "\t" + y) }
    pw.close
    println(timeStamp + "Constructing genome-architectures")
    //reload gene description and convert to list of paths
    val genome_paths = {
      descriptions.values.groupBy(_.ref).mapValues(_.toList.sortBy(x => (x.start, x.end))
        .map(x => new Gene(gene_families.getOrElse(x.id, x.id), x.ori)))
    }
    println(timeStamp + "--Constructed architectures for " + genome_paths.size + " sequences")
    println(timeStamp + "Constructing gene-graph with coverage statistics")
    //construct gene graph along with node and edge coverage
    val (graph, ncov, ecov) = {
      val tmp = genome_paths.foldLeft((empty_gene_graph, empty_node_coverage, empty_edge_coverage)) {
        case ((_graph, _ncov, _ecov), (seq, path)) => updateGeneGraph(path, false, _graph, _ncov, _ecov)
      }
      (tmp._1.mapValues(_.distinct), tmp._2, tmp._3)
    }
    println(timeStamp + "--Constructed gene graph with " + computeTotalNodes(graph) + " nodes and " + (ecov.size / 2) +
      " distinct edges")
    println(timeStamp + "Writing to disk")
    val pw2 = new PrintWriter(config.outputDir + "/" + config.prefix + ".gfa")
    ptolemyGraph2GFA(pw2, graph, genome_paths, ncov, ecov, Map(), Map())
    pw2.close
    println(timeStamp + "Successfully completed!")

  }

  def min: (Double, Double) => Double = (x, y) => if (x < y) x else y

}
